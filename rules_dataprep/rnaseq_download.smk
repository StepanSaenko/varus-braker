# Rule: retrieve_rnaseq_info_from_sra
# Purpose:
#   This rule automates the retrieval of RNA sequencing (RNA-seq) information for specified taxa
#   from the NCBI Sequence Read Archive (SRA). It utilizes a Python script to query NCBI for RNA-seq
#   data related to unannotated species listed in a provided table. The rule processes each taxon
#   listed in the 'data/checkpoints_dataprep/{taxon}_blank.tbl', extracting necessary information for downstream RNA-seq
#   analysis pipelines.
#
# Inputs:
#   - download_script: A Python script ('scripts/retrieve_rnaseq_info.py') that handles the querying
#     of NCBI databases and prepares output lists of SRA accession numbers and other relevant data.
#   - unannotated_species: A table listing species without annotated genomes; used to determine which
#     species' RNA-seq data to retrieve.
#
# Outputs:
#   - fastqdump_lst: A list file ('data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst') containing accession numbers
#     and other parameters formatted for the `prefetch` and `fasterq-dump` utility, which can be used to download sequence data.
#   - done: A simple checkpoint file ('data/checkpoints_dataprep/{taxon}_rnaseq_info.done') indicating successful completion of data retrieval.
#
# Parameters:
#   - taxon: A wildcard parameter defining the specific taxon being processed. Extracted from the filename pattern.
#   - email: An email address registered with NCBI, required for making API requests to NCBI databases.
#
# ToDo:
#   - make the number of libraries to download random & configurable in terms of number by config.ini
#   - make sure the libraries are mRNA and do not match microbiome or metagenome
rule retrieve_rnaseq_info_from_sra:
    input:
        download_script = "scripts/retrieve_rnaseq_info.py",
        unannotated_species = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl"
    output:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
    params:
        taxon = lambda wildcards: wildcards.taxon,
        email = config['ESEARCH']['email'],
        nlibs = int(config['FASTERQDUMP']['number_of_libraries'])
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        python3 {input.download_script} -e {params.email} -t {input.unannotated_species} -f {output.fastqdump_lst} -n {params.nlibs};
        """

# Rule: download_fastq
# Warning: 
#   This rule takes extremely long! It runs on SLURM but cannot take full advantage of a larger
#   number of threads since i/o is the bottleneck. We are using /dev/shm to speed up as much
#   possible.
#   It must regrettably expected that this rule fails, repeatedly. Do not worry, just restart the
#   workflow. The rule is implemented in a way that it will pick up where it left off. 
#   This rule is designed to process a list of SRA accession numbers and use `fasterq-dump` to download
#   paired-end fastq files for each accession. It handles the preparation of directory structures
#   and the organization of downloaded files within those directories. This is crucial for downstream
#   analysis such as RNA-seq data processing or genomic assemblies.
#
# Inputs:
#   - fastqdump_lst: A file ('data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst') that lists species and their
#     corresponding SRA accession numbers. Each line contains a species name and a comma-separated list
#     of SRA IDs.
#
# Outputs:
#   - done: A file ('data/checkpoints_dataprep/{taxon}_fastqdump.done') that signals the successful completion of the downloads
#     and processing of all SRA data listed in the input file.
#
# Operations:
#   - The rule first sets up the environment to ensure correct directory binding when using Singularity.
#   - It reads each line from the input file, processes the species name to replace spaces with underscores
#     (for file naming consistency), and creates a directory for storing the downloaded fastq files.
#   - For each SRA ID, `fasterq-dump` is executed to retrieve paired-end reads, and the resultant files are
#     immediately compressed using gzip to save space.
#   - A log file is generated to capture the output and errors of the download process, helping in troubleshooting
#     and ensuring transparency of the operation.
#   - Upon successful execution of all commands, a 'done' file is created as a checkpoint indicating
#     the completion of the task.
#
# This job uses at most 12 threads but occupies the entire RAM of a node because it stores data in RAM
# Rule: download_fastq
# Warning: 
#   This rule takes extremely long! It runs on SLURM but cannot take full advantage of a larger
#   number of threads since i/o is the bottleneck. We are using /dev/shm to speed up as much
#   possible.
#   It must regrettably expected that this rule fails, repeatedly. Do not worry, just restart the
#   workflow. The rule is implemented in a way that it will pick up where it left off. 
#   This rule is designed to process a list of SRA accession numbers and use `fasterq-dump` to download
#   paired-end fastq files for each accession. It handles the preparation of directory structures
#   and the organization of downloaded files within those directories. This is crucial for downstream
#   analysis such as RNA-seq data processing or genomic assemblies.
#
# Inputs:
#   - fastqdump_lst: A file ('data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst') that lists species and their
#     corresponding SRA accession numbers. Each line contains a species name and a comma-separated list
#     of SRA IDs.
#
# Outputs:
#   - done: A file ('data/checkpoints_dataprep/{taxon}_fastqdump.done') that signals the successful completion of the downloads
#     and processing of all SRA data listed in the input file.
#
# Operations:
#   - The rule first sets up the environment to ensure correct directory binding when using Singularity.
#   - It reads each line from the input file, processes the species name to replace spaces with underscores
#     (for file naming consistency), and creates a directory for storing the downloaded fastq files.
#   - For each SRA ID, `fasterq-dump` is executed to retrieve paired-end reads, and the resultant files are
#     immediately compressed using gzip to save space.
#   - A log file is generated to capture the output and errors of the download process, helping in troubleshooting
#     and ensuring transparency of the operation.
#   - Upon successful execution of all commands, a 'done' file is created as a checkpoint indicating
#     the completion of the task.
#
# This job uses at most 12 threads but occupies the entire RAM of a node because it stores data in RAM
rule run_download_fastq:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_A09_shorten_genomic_headers.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B03_fastqdump.done"
    params:
        threads = 10,
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    threads: 12
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        echo "I am done with the SIF file, now starting, will take very long..."
        logfile=$PWD/data/checkpoints_dataprep/{params.taxon}_B03_fastqdump.log
        # Read the input file and process it
        while IFS=$'\\t' read -r species sra_ids; do
            echo "$PWD" &>> $logfile
            species_fixed=$(echo "$species" | sed 's/ /_/g')  # Replace space with underscore
            # if the directory fastq does not exist, yet
            if [ ! -d "data/species/$species_fixed/fastq" ]; then
                echo "mkdir -p data/species/$species_fixed/fastq" &>> $logfile
                mkdir -p data/species/$species_fixed/fastq  # Create a directory for the species
            fi
            if [ ! -d "data/species/$species_fixed/sra" ]; then
                echo "mkdir -p data/species/$species_fixed/sra" &>> $logfile
                mkdir -p data/species/$species_fixed/sra  # Create a directory for the species
            fi
            # Convert comma-separated string to array
            IFS=',' read -ra ids <<< "$sra_ids"

            # Process each SRA ID locally using the array
            for id in "${{ids[@]}}"; do
                # this is such a long and expensive process that we do not want it to execute if fastq.gz files already exist
                if [ ! -f "data/species/$species_fixed/fastq/${{id}}_1.fastq.gz" ] && [ ! -f "data/species/$species_fixed/fastq/${{id}}_2.fastq.gz" ]; then
                    echo "prefetch -O data/species/$species_fixed/sra $id" &>> $logfile
                    prefetch -O data/species/$species_fixed/sra $id --max-size 30g &>> $logfile
                    echo "fasterq-dump data/species/$species_fixed/sra/$id/$id.sra --split-files -O data/species/$species_fixed/fastq -e {params.threads}" &>> $logfile
                    fasterq-dump data/species/$species_fixed/sra/$id/$id.sra --split-files -O data/species/$species_fixed/fastq -e {params.threads} &>> $logfile
                    echo "gzip data/species/$species_fixed/fastq/${{id}}_1.fastq" &>> $logfile
                    gzip data/species/$species_fixed/fastq/${{id}}_1.fastq &>> $logfile
                    gzip data/species/$species_fixed/fastq/${{id}}_2.fastq
                fi
            done
            rm -rf data/species/$species_fixed/sra
        done < {input.fastqdump_lst} &>> $logfile
        
        echo "touch {output.done}" &>> $logfile
        touch {output.done}
        """


# Rule 'run_hisat2_index' prepares a HISAT2 index for RNA-Seq data by leveraging the hisat2-build command. 
# It ensures that all required inputs are available including the RNA-Seq data list, signal files indicating 
# successful download and genome preparation. Outputs a completion flag upon successfully building the index.
rule run_hisat2_index:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        download_done = "data/checkpoints_dataprep/{taxon}_B03_fastqdump.done",
        genome_done = "data/checkpoints_dataprep/{taxon}_A09_shorten_genomic_headers.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B04_hisat2_index.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B04_hisat2_index.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            if [ ! -d "data/species/$species/hisat2" ]; then
                mkdir -p "data/species/$species/hisat2"
            fi
            echo "hisat2-build -p {params.threads} data/species/$species/genome/genome.fa data/species/$species/genome/genome" >> $log
            which hisat2-build &>> $log
            hisat2-build -p {params.threads} data/species/$species/genome/genome.fa data/species/$species/genome/genome &>> $log
        done
        touch {output.done}
        """


# Rule 'run_hisat2' performs RNA sequencing alignment using HISAT2. It processes each RNA-Seq dataset listed
# in the fastqdump list, provided the HISAT2 indexing has been completed. The rule generates alignment files (.sam) 
# for each dataset and logs the process details. It operates within a specified Singularity container to ensure
# compatibility and consistency of the software environment.
rule run_hisat2:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B04_hisat2_index.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B05_hisat2.done",
        log = "data/checkpoints_dataprep/{taxon}_B05_hisat2.log"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        log={output.log}
        # if file $log does not exist yet
        if [ ! -f $log ]; then
            touch $log
        fi
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.sam" ]; then
                    echo "hisat2 -p {params.threads} -x data/species/$species/genome/genome -1 data/species/$species/fastq/${{sra_id}}_1.fastq.gz -2 data/species/$species/fastq/${{sra_id}}_2.fastq.gz -S data/species/$species/hisat2/${{sra_id}}.sam" >> $log
                    hisat2 -p {params.threads} -x data/species/$species/genome/genome -1 data/species/$species/fastq/${{sra_id}}_1.fastq.gz -2 data/species/$species/fastq/${{sra_id}}_2.fastq.gz -S data/species/$species/hisat2/${{sra_id}}.sam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sam already exists" &>> $log
                fi
            done
        done
        touch {output.done}
        """

# During identification of RNA-Seq libraries from SRA, we do not exclude co-culture libraries
# because in taxa, such as diatoms, most data are obtained from such libraries. However, we
# thereby also download libraries that are not mRNA libraries. To exclude these libraries, we
# use the hisat2 log file to identify libraries with less than a certain percentage of reads
# that are aligned. We then remove these libraries from the list of libraries to process.
# ToDo: This possibly leaves dangling large sam file that need to be removed, later!
# They are not removed, here, to avoid repeated download & alignment in case the pipeline needs 
# to be relaunched.
rule remove_bad_libraries:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        hisat2_log = "data/checkpoints_dataprep/{taxon}_B05_hisat2.log"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B06_remove_bad_libraries.done",
        new_fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
    params:
        taxon=lambda wildcards: wildcards.taxon,
        mapping_threshold = float(config['FASTERQDUMP']['minimial_aligned_reads_percent']),
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        '''
        log=data/checkpoints_dataprep/{params.taxon}_B06_remove_bad_libraries.log
        cmd="python3 scripts/remove_bad_sra.py -t {params.mapping_threshold} -l {input.hisat2_log} -o {output.new_fastqdump_lst} -p data/species -i {input.fastqdump_lst}"
        echo $cmd &> $log
        $cmd &>> $log
        touch {output.done}
        '''


# Rule 'run_sam_to_bam' converts SAM files to BAM format using Samtools, optimizing for reduced file size and
# enhanced processing speeds. It only processes samples that have completed the 'remove_bad_libraries' step, ensuring 
# only quality data is transformed. The rule generates a log of the conversion process and a signal file upon completion.
rule run_sam_to_bam:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        remove_bad_done = "data/checkpoints_dataprep/{taxon}_B06_remove_bad_libraries.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B06_sam2bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B06_sam2bam.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.bam" ] && [ -f data/species/$species/hisat2/${{sra_id}}.sam ]; then
                    echo "samtools view --threads {params.threads} -bS data/species/$species/hisat2/${{sra_id}}.sam > data/species/$species/hisat2/${{sra_id}}.bam" &>> $log
                    samtools view --threads {params.threads} -bS data/species/$species/hisat2/${{sra_id}}.sam > data/species/$species/hisat2/${{sra_id}}.bam 2>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.bam already exists" &>> $log
                fi
            done
        done
        touch {output.done}
        """


# Rule 'run_samtools_sort_single' sorts BAM files by genomic coordinates using Samtools, 
# which facilitates more efficient downstream analyses such as feature counting. This 
# rule acts upon BAM files generated in prior steps once they are confirmed complete 
# by the 'sam2bam.done' checkpoint. It logs the sorting process and produces a sorted 
# BAM for each RNA-Seq sample, enhancing the organization and accessibility of sequence 
#data.
rule run_samtools_sort_single:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B06_sam2bam.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B07_samtools_sort_single.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B07_samtools_sort_single.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.sorted.bam" ] && [ -f data/species/$species/hisat2/${{sra_id}}.bam ]; then
                    echo "samtools sort --threads {params.threads} data/species/$species/hisat2/${{sra_id}}.bam -o data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                    samtools sort --threads {params.threads} data/species/$species/hisat2/${{sra_id}}.bam -o data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sorted.bam already exists" &>> $log
                fi
            done
        done
        touch {output.done}
        """


# Rule 'run_samtools_index_single' indexes sorted BAM files using Samtools, creating 
# an index file for each. This indexing step is crucial for efficient retrieval of 
# data during downstream genomic analyses such as read depth calculations or localized 
# querying. The rule triggers only after the sorting of BAM files is confirmed complete, 
# ensuring a streamlined workflow. Outputs include a completion signal for the indexing 
# process.
rule run_samtools_index_single:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B07_samtools_sort_single.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B08_samtools_index_single.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B08_samtools_index_single.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                echo "samtools index data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                samtools index data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
            done
        done
        touch {output.done}
        """


# Rule 'cleanup_sam_bam_unsorted_files' is responsible for cleaning up intermediate 
# SAM and unsorted BAM files after indexing is completed to conserve disk space. 
# This rule is triggered only after the successful completion of BAM indexing, 
# ensuring that only redundant files are removed. The process is logged for audit
# purposes, and a completion file is generated once all specified deletions are executed. 
# This cleanup step helps maintain an organized file system and reduces storage overhead,
# essential in large-scale sequencing projects.
rule cleanup_sam_bam_unsorted_files:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B08_samtools_index_single.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B09_cleanup_sam_bam_unsorted.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B09_cleanup_sam_bam_unsorted.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            echo "Processing species $species" &>> $log
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                echo "Processing SRA ID $sra_id" &>> $log
                if [ -f data/species/$species/hisat2/${{sra_id}}.sam ]; then
                    echo "rm data/species/$species/hisat2/${{sra_id}}.sam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.sam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sam does not exist" &>> $log
                fi
                if [ -f data/species/$species/hisat2/${{sra_id}}.bam ]; then
                    echo "rm data/species/$species/hisat2/${{sra_id}}.bam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.bam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.bam does not exist" &>> $log
                fi
            done
        done
        touch {output.done}
        """


# Rule 'run_merge_bam' is designed to consolidate sorted BAM files for each species 
# into a single BAM file, facilitating streamlined downstream analyses. This 
# merging process is triggered only after all unsorted and intermediate BAM and 
# SAM files have been cleaned up. The rule logs detailed processing steps and 
# conditions, such as the number of files merged and any copying actions for 
# species with a single sorted BAM file. It utilizes multi-threading to optimize 
# the merging process.
rule run_merge_bam:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B09_cleanup_sam_bam_unsorted.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B10_merge_bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B10_samtools_merge.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            echo "Processing species $species" &>> $log
            # count number of files matching the pattern data/species/$species/hisat2/*.sorted.bam, excluding file data/species/$species/hisat2/$species.sorted.bam
            num_files=$(ls -1 data/species/$species/hisat2/*.sorted.bam | grep -v data/species/$species/hisat2/$species.sorted.bam | wc -l)
            echo "Number of files is $num_files" &>> $log
            # merge all bam files of the same species with samtools merge 
            if [ $num_files -gt 1 ] && [ ! -f "data/species/$species/hisat2/$species.bam" ]; then
                echo "samtools merge --threads {params.threads} -o data/species/$species/hisat2/$species.sorted.bam data/species/$species/hisat2/*.sorted.bam" &>> $log
                samtools merge --threads {params.threads} -o data/species/$species/hisat2/$species.bam data/species/$species/hisat2/*.sorted.bam &>> $log
            elif [ $num_files -eq 1 ] && [ ! -f "data/species/$species/hisat2/$species.bam" ] ; then
                # get the actual file name from wild card expansion
                filename=$(echo data/species/$species/hisat2/*.sorted.bam)
                echo "mv $filename data/species/$species/hisat2/$species.bam" &>> $log
                mv $filename data/species/$species/hisat2/$species.bam &>> log
            else
                echo "data/species/$species/hisat2/$species.bam already exists" &>> $log
            fi
        done
        touch {output.done}
        """


# Rule 'cleanup_sorted_bam_files' handles the removal of individual sorted BAM files
# once they have been successfully merged into a single BAM file per species. This 
# step is critical for managing storage efficiently, especially in large-scale genomic
# projects where data volume can become unwieldy. The rule executes only after the 
# merging process is verified complete, ensuring that no data necessary for further
# analysis is inadvertently lost. It logs each deletion step for traceability and 
# issues a completion signal once all specified files are removed. This cleanup is 
# essential for maintaining an organized data structure and minimizing disk space usage.
rule cleanup_sorted_bam_files:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B10_merge_bam.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B11_cleanup_sorted_bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B11_cleanup_single_bam.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            # delete the individual bam files
            for sra_id in "${{sra_array[@]}}"; do
                if [ -f "data/species/$species/hisat2/${{sra_id}}.sorted.bam" ]; then
                    echo "rm data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sorted.bam does not exist" &>> $log
                fi
            done
        done
        touch {output.done}
        """

# Rule 'run_sort_merged_bam' further processes the merged BAM files by sorting 
# them again using Samtools to ensure that they are optimized for downstream 
# genomic analyses. This step is crucial for the efficient handling of large BAM 
# files, particularly in preparation for tasks such as variant calling or 
# statistical analyses. The rule activates only after the successful completion 
# of the BAM merging process, guaranteeing that all necessary preliminary 
# steps have been fulfilled. It employs significant computational resources 
# to handle the demands of sorting large genomic datasets and logs each operation 
# for tracking and debugging purposes.
rule run_sort_merged_bam:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        merge_done = "data/checkpoints_dataprep/{taxon}_B10_merge_bam.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B12_sort_merged_bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B12_samtools_sort_merged.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            # merge all bam files of the same species with samtools merge 
            if [ ! -f "data/species/$species/hisat2/${{species}}.sorted.bam" ] && [ -f  data/species/$species/hisat2/${{species}}.bam ]; then
                echo "samtools sort --threads {params.threads} data/species/$species/hisat2/${{species}}.bam -o data/species/$species/hisat2/${{species}}.sorted.bam &>> $log" &>> $log
                samtools sort --threads {params.threads} data/species/$species/hisat2/${{species}}.bam -o data/species/$species/hisat2/rnaseq.s.bam &>> $log
            else
                echo "data/species/$species/hisat2/${{species}}.sorted.bam already exists" &>> $log
            fi
        done
        touch {output.done}
        """

rule cleanup_rnaseq:
    input:
        sort_merged = "data/checkpoints_dataprep/{taxon}_B12_sort_merged_bam.done",
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B13_cleanup_rnaseq.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        log="data/checkpoints_dataprep/{params.taxon}_B13_cleanup.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            if [ -d "data/species/$species/fastq"]; then
                echo "rm -r data/species/$species/fastq" &>> $log
                rm -r data/species/$species/fastq &>> $log
            fi
            # find all sam files in data/species/$species/hisat2 and store them in an array
            readarray -t sam_files < <(find data/species/$species/hisat2 -name "*.sam")
            for sam_file in "${{sam_files[@]}}"; do
                echo "rm $sam_file" &>> $log
                rm $sam_file &>> $log
            done
            # implement removal of hisat2 index!
        done
        touch {output.done}
        """

'''
# A perfectly fine VARUS rule that does not work simply because VARUS tried to chdir
# chdir is not allowed in snakemake.
rule run_varus:
    input:
        varus_lst = "data/checkpoints_dataprep/{taxon}_rnaseq_for_varus.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_download.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_varus.done"
    params:
        threads = int(config['SLURM_ARGS']['cpus_per_task']),
        taxon=lambda wildcards: wildcards.taxon
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    threads: 12
    resources:
        mem_mb=int(int(config['SLURM_ARGS']['mem_of_node'])/int(config['SLURM_ARGS']['cpus_per_task'])*12),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        log="data/checkpoints_dataprep/{params.taxon}_varus.log"
        echo "" > $log
        readarray -t lines < <(cat {input.varus_lst})
        for line in "${{lines[@]}}"; do
            echo "Original line: $line" >> $log
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            echo "Modified species line: $species" >> $log
            pwd >> $log
            if [ ! -d "data/species/$species/varus" ]; then
                mkdir -p "data/species/$species/varus"
            fi
            IFS='_' read -r genus spec <<< "$species"
            echo "runVARUS.pl --aligner=HISAT --runThreadN={params.threads} --speciesGenome=../genome/genome.fa --readFromTable=0 --createindex=1 --verbosity=5 --latinGenus=$genus --latinSpecies=$spec --varusParameters=VARUSparameters.txt --outFileDir data/species/$species/varus 2> data/species/$species/varus/varus.err" &>> $log
            ( runVARUS.pl --aligner=HISAT --runThreadN={params.threads} \
                --speciesGenome=../genome/genome.fa \
                --readFromTable=0 --createindex=1 --verbosity=5 \
                --latinGenus=$genus --latinSpecies=$spec \
                --varusParameters=VARUSparameters.txt \
                --outFileDir data/species/$species/varus \
                2> data/species/$species/varus/varus.err )
        done
        touch {output.done}
        """
'''
