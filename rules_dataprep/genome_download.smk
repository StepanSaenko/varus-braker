# Rule: download_assembly_info
# Purpose:
# This rule automates the process of downloading genome assembly info/ for a specified taxon
# from the NCBI database. It utilizes the NCBI datasets command-line tool within a Singularity
# container to ensure environment consistency and reproducibility. The downloaded data is
# initially compressed into a ZIP file, which is then extracted, and the relevant JSONL file
# containing assembly data is moved to a specified directory.
#
# Inputs:
#   - None directly specified, but the rule uses a parameter 'taxon' derived from wildcards,
#     typically provided by upstream rules or workflow orchestration.
#
# Outputs:
#   - raw_json: A JSONL file containing detailed assembly information for the specified taxon.
#     This file is used in downstream analysis and must be stored in the 'data' directory.
#
# Singularity:
#   - Uses a Docker container 'katharinahoff/varus-notebook' converted to Singularity
#     image format, ensuring that the NCBI datasets tool and all its dependencies are correctly
#     configured and isolated from the host environment.
#
# Steps Executed:
# 1. Set up a bind mount to ensure current working directory ('${PWD}') is accessible inside the Singularity container.
# 2. Execute the 'datasets download genome' command specifying the taxon, source database (GenBank), and download format (dehydrated).
# 3. Unzip the downloaded file to a specified sub-dirlocalrules: all, another_local_ruleectory.
# 4. Move the assembly data report (JSONL format) to the designated output location.
# 5. Clean up all intermediate files and directories to maintain a clean working environment.
rule download_assembly_info:
    output:
        raw_json = "data/checkpoints_dataprep/{taxon}_A01.json"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        datasets download genome taxon "{params.taxon}" --assembly-source genbank --dehydrated --filename {params.taxon}_ncbi.zip; \
        unzip -o {params.taxon}_ncbi.zip -d {params.taxon}_ncbi_dataset; \
        mv {params.taxon}_ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl {output.raw_json}; \
        rm -rf {params.taxon}_ncbi_dataset {params.taxon}_ncbi.zip; \
        mkdir -p data/species
        """


# Rule: assembly_json_to_tbl
# Purpose:
# This rule processes JSON files containing genome assembly metadata for different taxa.
# Each JSON file may contain one or multiple JSON objects representing different assembly reports.
# The rule reads the JSON file, parses multiple JSON objects if present, and then extracts
# relevant fields to produce a tab-separated values (TSV) file. The output file includes columns for
# accession, species name, assembly status, protein coding gene counts, contig N50 size, and the refSeq category.
# This rule handles JSON reading and parsing errors and will raise exceptions if files cannot be read or written.
#
# Input:
#   json_file - Path to the JSON file for a specific taxon, located in the 'data' directory.
# Output:
#   processed_tbl - Path to the output TSV file which will also be located in the 'data' directory.
# Processing Steps:
# 1. Read the entire JSON file content.
# 2. Split the content into individual JSON strings based on a specified delimiter.
# 3. Convert each JSON string back into a JSON object (dictionary) and handle any JSON parsing errors.
# 4. Extract the required information from each JSON object and write it to the output TSV file.
# 5. Handle and report any file input/output errors during reading and writing operations.
rule assembly_json_to_tbl:
    input:
        json_file = "data/checkpoints_dataprep/{taxon}_A01.json"
    output:
        processed_tbl = "data/checkpoints_dataprep/{taxon}_A02.tbl"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        taxon = wildcards.taxon
        print(taxon)
        json_file_path = f"data/checkpoints_dataprep/{taxon}_A01.json"
        tbl_file_path = f"data/checkpoints_dataprep/{taxon}_A02.tbl"
        try:
            print(json_file_path)
            # Read the entire file into a single string
            with open(json_file_path, 'r') as file:
                file_content = file.read()
        except IOError:
            raise Exception(f"File not found: {json_file_path}")
        # Split the content by a delimiter that indicates the end of a JSON object, assuming each JSON object ends with '}' on a new line
        json_objects = file_content.split('}\n')
        # Filter out empty strings if any and add closing brace which was removed by split
        json_objects = [obj + '}' for obj in json_objects if obj.strip()]
        data = []
        for json_str in json_objects:
            try:
                # Load each JSON object as a dictionary
                data.append(json.loads(json_str))
            except json.JSONDecodeError as e:
                print(f"Failed to decode JSON: {e}")
            continue
        try:
            with open(tbl_file_path, "w") as f:
                f.write("accession\tspecies\tstatus\tproteins\tcontigN50\trefseqCategory\n")
                # If data is a list of entries, iterate through it
                # If it's a single dict, adjust the loop or wrap `data` in a list: [data]
                for entry in data:
                    accession = entry.get('accession', '')
                    species = entry['organism'].get('organismName', '')
                    # only proceed if species does not match sp. as part of the string
                    if "sp." in species or "uncultured" in species:
                        continue
                    if len(species.split(" ")) != 2:
                        # cut off everthing after the first two words
                        species = " ".join(species.split(" ")[:2])
                    status = entry['assemblyInfo'].get('assemblyStatus', '')
                    proteins = entry.get('annotationInfo', {}).get('stats', {}).get('geneCounts', {}).get('proteinCoding', 'N/A')
                    contigN50 = entry['assemblyStats'].get('contigN50', '')
                    refSeqCategory = entry['assemblyInfo'].get('refseqCategory', 'N/A')
                    f.write(f"{accession}\t{species}\t{status}\t{proteins}\t{contigN50}\t{refSeqCategory}\n")
        except IOError:
            raise Exception(f"Error writing to file: {tbl_file_path}")


# Rule: classify_species
# Purpose:
# This rule is designed to classify species based on their protein counts from genome assembly data.
# It separates the species into two categories: those that are well-annotated (with protein counts
# exceeding 1000) and those that are either poorly annotated or have missing protein data ('N/A').
#
# Inputs:
#   - tbl_file: A tab-separated file ('.tbl') for each taxon, containing genome assembly data.
#     This file includes various metrics among which are protein counts.
#
# Outputs:
#   - already_annotated_tbl: Outputs a file listing species that are considered well-annotated.
#   - not_annotated_tbl: Outputs a file listing species that are either not well-annotated or
#     whose protein count data is missing.
#
# Parameters:
#   - taxon: A wildcard parameter that dynamically accepts the name of the taxon from the workflow.
#
# Processing Steps:
# 1. Read the input TBL file into a pandas DataFrame. Non-numeric 'proteins' values are converted to NaN
#    to facilitate numerical operations and comparisons.
# 2. Filter the DataFrame to create two subsets:
#    - 'annotated_data': Contains entries with protein counts greater than 1000.
#      This subset is further processed to keep only the highest quality genome per species,
#      prioritized by 'refseqCategory' and the highest protein count.
#    - 'blank_data': Contains entries with protein counts less than or equal to 1000 or with missing protein data.
#      This subset is sorted by 'contigN50' to prioritize genomes based on genome assembly quality, keeping the
#      highest 'contigN50' value per species.
# 3. Both subsets are written to separate output files in a tab-separated format.
rule classify_species:
    input:
        tbl_file = "data/checkpoints_dataprep/{taxon}_A02.tbl"
    output:
        already_annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl",
        not_annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        taxon = wildcards.taxon
        tbl_file_path = f"data/checkpoints_dataprep/{taxon}_A02.tbl"
        annotated_tbl_path = f"data/checkpoints_dataprep/{taxon}_A03_annotated.tbl"
        blank_tbl_path = f"data/checkpoints_dataprep/{taxon}_A03_blank.tbl"
        
        # Read the data, handling non-numeric proteins values
        try:
            data = pd.read_csv(tbl_file_path, sep="\t")
            data['proteins'] = pd.to_numeric(data['proteins'], errors='coerce')  # Convert 'N/A' to NaN
        except IOError:
            raise Exception(f"Error reading file: {tbl_file_path}")

        # Make a subset of representative genomes, this can be seen by refseqCategory being not empty
        representative_genomes = data[data['refseqCategory'].notna()]
        # Sort this subset by the highest protein number
        representative_genomes = representative_genomes.sort_values(by=['proteins'], ascending=False).drop_duplicates(subset='species', keep='first')
        # go through representative genomes and check whether there are species duplicates, if yes, only keep the first
        representative_genomes = representative_genomes.drop_duplicates(subset='species', keep='first')
        # from data, get only those that are na in the refseq Cateogry, store in non_representative_genomes
        non_representative_genomes = data[data['refseqCategory'].isna()]
        # sort the non_representative_genomes by protein number
        non_representative_genomes = non_representative_genomes.sort_values(by=['proteins'], ascending=False).drop_duplicates(subset='species', keep='first')
        data = representative_genomes
        # if data has no rows, create an empty data frame data that has exacly the same columns 
        if data.empty:
            data = pd.DataFrame(columns=['accession', 'species', 'status', 'proteins', 'contigN50', 'refseqCategory'])

        # add to all_filtered genomes those lines from non_represenative genomes that have no species entry in all_filtered_genomes yet
        for index, row in non_representative_genomes.iterrows():
            if row['species'] not in data['species'].values:
                data = pd.concat([data, pd.DataFrame([row])], ignore_index=True)

        # Annotated data filter
        annotated_data = data[data['proteins'] > 1000]
        
        # Blank data filter
        blank_data = data[(data['proteins'] <= 1000) | (data['proteins'].isna())]

        # Output data to files
        try:
            annotated_data.to_csv(annotated_tbl_path, sep="\t", index=False)
            blank_data.to_csv(blank_tbl_path, sep="\t", index=False)
        except IOError:
            raise Exception(f"Error writing to file: {annotated_tbl_path} or {blank_tbl_path}")


# This Snakemake rule 'prepare_download_assemblies_from_ncbi' is designed to automate the preparation of shell scripts
# for downloading genomic assemblies from the NCBI database. The rule uses specific table files containing species and 
# accession numbers to generate commands for downloading and organizing genomic data into structured directories. The
# output is a shell script tailored for each taxon that handles directory creation, data retrieval, and data extraction.
rule prepare_download_assemblies_from_ncbi:
    input:
        annotated_tbl_file = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl",
        blank_tbl_file = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl"
    output:
        download_script = "data/checkpoints_dataprep/{taxon}_A04_download.sh"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        # Read the annotated table to get accessions and species names
        print(input.blank_tbl_file)
        df_blank = pd.read_csv(input.blank_tbl_file, sep="\t", usecols=['accession', 'species'])
        df_anno = pd.read_csv(input.annotated_tbl_file, sep="\t", usecols=['accession', 'species'])
        # Prepare a file to capture assembly information
        try:
            with open(output.download_script, 'w') as outfile:
                command = ""
                for index, row in df_blank.iterrows():
                    # take the row["species"] field and replace all spaces by _, otherwise remove all special characters
                    species = row["species"].replace(" ", "_")
                    species = re.sub(r'[/?,.*&\\;]+', '', species)
                    command += f"cd data/species; "
                    command += f"mkdir {species}; cd {species};"
                    command += f"datasets download genome accession {row['accession']} --filename {row['accession']}_assembly.zip; "
                    command += f"unzip -o {row['accession']}_assembly.zip; "
                    command += f"mkdir genome; "
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.fna genome/genome.fa; "
                    command += f"rm -rf ncbi_dataset {row['accession']}_assembly.zip; cd ../../..;\n"
                for index, row in df_anno.iterrows():
                    # take the row["species"] field and replace all spaces by _, otherwise remove all special characters
                    species = row["species"].replace(" ", "_")
                    species = re.sub(r'[/?,.*&\\;]+', '', species)
                    command += f"cd data/species; "
                    command += f"mkdir {species}; cd {species};"
                    command += f"datasets download genome accession {row['accession']} --filename {row['accession']}_assembly.zip --include genome,gff3,protein,gtf; "
                    command += f"unzip -o {row['accession']}_assembly.zip; "
                    command += f"mkdir genome; "
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.fna genome/genome.fa; "
                    command += f"mkdir annot;"
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.gff annot/annot.gff3;"
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.gtf annot/annot.gtf;"
                    command += f"mkdir prot;"
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.faa prot/protein.faa;"
                    command += f"rm -rf ncbi_dataset {row['accession']}_assembly.zip; rm README.md; cd ../../..;"
                outfile.write(command)
        except IOError:
            raise Exception(f"Error writing to file: {output.download_script}")


# This rule executes the download commands specified in the download script generated by the rule 'prepare_download_assemblies_from_ncbi'. 
# The rule is not parallelized as the data input/output operations are likely to be the bottleneck.
rule execute_genome_download_commands:
    input:
        download_script = "data/checkpoints_dataprep/{taxon}_A04_download.sh"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A04_download.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        bash {input.download_script}; \
        touch {output.done}
        """

# If a different assembly of a species was already annotated, and if we have these protein
# sequences, we may use them to annotate the new assembly. This rule prepares a shell script
# that downloads the protein sequences of the already annotated species from the NCBI database.
rule prepare_legacy_protein_download:
    input:
        all_tbl = "data/checkpoints_dataprep/{taxon}_A02.tbl",
        not_annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl",
        base_download_done = "data/checkpoints_dataprep/{taxon}_A04_download.done"
    output:
        legacy_proteins_script = "data/checkpoints_dataprep/{taxon}_A05_legacy_proteins_prep.sh"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        df_all = pd.read_csv(input.all_tbl, sep="\t", usecols=['accession', 'species', 'proteins', 'contigN50'])
        # drop all rows that do not have a number >1000 in proteins column
        # this drops the organelle only annotations
        df_all = df_all[df_all['proteins'] > 1000]
        # if duplicates in species exist, dort by contigN50 and keep only the rows with the longest contigN50
        df_all = df_all.sort_values(by=['contigN50'], ascending=False).drop_duplicates(subset='species', keep='first')
        df_blank = pd.read_csv(input.not_annotated_tbl, sep="\t", usecols=['accession', 'species'])
        # for every row in df_blank check whether the species exists in df_all, and whether the accession number is different
        # if yes, print the accession number from df_all
        try:
            with open(f"data/checkpoints_dataprep/{params.taxon}_A05_legacy_proteins_prep.sh", 'w') as outfile:
                command = ""
                for index, row in df_blank.iterrows():
                    if row['species'] in df_all['species'].values:
                        if row['accession'] != df_all[df_all['species'] == row['species']]['accession'].values[0]:
                            command += f"if [ ! -f data/species/{row['species'].replace(' ', '_')}/prot_legacy/proteins.faa ] && [ ! -f data/species/{row['species'].replace(' ', '_')}/prot_legacy/proteins.fa ]; then\n"
                            command += f"\tmkdir data/species/{row['species'].replace(' ', '_')}/prot_legacy;\n"
                            command += f"\tcd data/species/{row['species'].replace(' ', '_')}/prot_legacy;\n"
                            command += f"\tdatasets download genome accession {df_all[df_all['species'] == row['species']]['accession'].values[0]} --filename {row['accession']}_legacy_proteins.zip --include protein;\n"
                            command += f"\tunzip -o {row['accession']}_legacy_proteins.zip;\n"
                            command += f"\tmv ncbi_dataset/data/{df_all[df_all['species'] == row['species']]['accession'].values[0]}/*.faa proteins.faa;\n"
                            command += f"\trm -rf ncbi_dataset {row['accession']}_legacy_proteins.zip; rm README.md; cd ../../../..;\n"
                            command += f"fi\n"
                outfile.write(command)
        except IOError:
            raise Exception(f"Error writing to file: data/checkpoints_dataprep/{wildcards.taxon}_A05_legacy_proteins_prep.sh")


rule execute_legacy_prot_download_commands:
    input:
        download_script = "data/checkpoints_dataprep/{taxon}_A05_legacy_proteins_prep.sh"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A06_legacy_proteins_download.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        bash {input.download_script}; \
        touch {output.done}
        """

rule write_simplify_legacy_protein_headers_cmds:
    input:
        legacy_present = "data/checkpoints_dataprep/{taxon}_A06_legacy_proteins_download.done",
        all_tbl = "data/checkpoints_dataprep/{taxon}_A02.tbl",
        not_annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl"
    output:
        simplified_legacy_proteins = "data/checkpoints_dataprep/{taxon}_A07_simplify_legacy_proteins.sh"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        df_all = pd.read_csv(input.all_tbl, sep="\t", usecols=['accession', 'species', 'proteins', 'contigN50'])
        # drop all rows that do not have a number >1000 in proteins column
        # this drops the organelle only annotations
        df_all = df_all[df_all['proteins'] > 1000]
        # if duplicates in species exist, dort by contigN50 and keep only the rows with the longest contigN50
        df_all = df_all.sort_values(by=['contigN50'], ascending=False).drop_duplicates(subset='species', keep='first')

        df_blank = pd.read_csv(input.not_annotated_tbl, sep="\t", usecols=['accession', 'species'])
        # for every row in df_blank check whether the species exists in df_all, and whether the accession number is different
        # if yes, print the accession number from df_all
        try:
            with open(f"data/checkpoints_dataprep/{params.taxon}_A07_simplify_legacy_proteins.sh", 'w') as outfile:
                command = ""
                for index, row in df_blank.iterrows():
                    if row['species'] in df_all['species'].values:
                        if row['accession'] != df_all[df_all['species'] == row['species']]['accession'].values[0]:
                            command += f"if [ ! -f data/species/{row['species'].replace(' ', '_')}/prot_legacy/proteins.fa ]; then\n"
                            command += f"\tcd data/species/{row['species'].replace(' ', '_')}/prot_legacy; "
                            command += f"\tsimplifyFastaHeaders.pl proteins.faa prot_ proteins.fa header.map;"
                            command += f"\trm proteins.faa; cd ../../../..;\n"
                            command += f"fi\n"
                outfile.write(command)
        except IOError:
            raise Exception(f"Error writing to file: data/checkpoints_dataprep/{wildcards.taxon}_A06_legacy_prot_headers.sh")


rule execute_fix_headers_commands:
    input:
        download_script = "data/checkpoints_dataprep/{taxon}_A07_simplify_legacy_proteins.sh"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A08_fixed_protein_headers.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        bash {input.download_script}; \
        touch {output.done}
        """

# shorten both the genomic and the protein fasta headers (annotated proteins only, not the legacy proteins)
rule shorten_genomic_fasta_headers:
    input:
        genomic_download = "data/checkpoints_dataprep/{taxon}_A04_download.done",
        annotated_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl",
        blank_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl",
        legacy_proteins_download = "data/checkpoints_dataprep/{taxon}_A08_fixed_protein_headers.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A09_shorten_genomic_headers.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        df_anno = pd.read_csv(input.annotated_tbl_path, sep="\t", usecols=['species']) # here also proteins
        df_blank = pd.read_csv(input.blank_tbl_path, sep="\t", usecols=['species'])
        # concatenate the species of df_anno to df_blank
        df_blank = pd.concat([df_anno, df_blank]) # here genome
        # build paths to genome files from df_anno species field, where the space must be replaced by underscore
        genome_files = [f"data/species/{species.replace(' ', '_')}/genome/genome.fa" for species in df_blank['species']]
        # move all these files to genome.fa_original
        for genome_file in genome_files:
            shutil.move(genome_file, genome_file + "_original")
            # store the new file name as string
            original_genome_file = genome_file + "_original"
            try:
                with open(original_genome_file, 'r') as original_genome_handle:
                    # if a line starts with ">" remove everthing after the first whitespace, including the whitespace
                    # write the new line to a new file
                    with open(genome_file, 'w') as new_genome_handle:
                        for line in original_genome_handle:
                            if line.startswith(">"):
                                new_genome_handle.write(line.split(" ")[0] + "\n")
                            else:
                                new_genome_handle.write(line)
            except IOError:
                raise Exception(f"Error reading or writing to file: {original_genome_file} or {genome_file}")
            # delete the original file
            os.remove(original_genome_file)

            
        # build paths to protein files from df_anno species field, where the space must be replaced by underscore
        protein_files = [f"data/species/{species.replace(' ', '_')}/prot/protein.faa" for species in df_anno['species']]
        # move all these files to protein.faa_original
        for protein_file in protein_files:
            shutil.move(protein_file, protein_file + "_original")
            # store the new file name as string
            original_protein_file = protein_file + "_original"
            try:
                with open(original_protein_file, 'r') as protein_handle:
                    # if a line starts with ">" remove everthing after the first whitespace, including the whitespace, replace possible dots by underscors
                    # write the new line to a new file
                    with open(protein_file, 'w') as new_protein_handle:
                        for line in protein_handle:
                            if line.startswith(">"):
                                new_protein_handle.write(line.split(" ")[0].replace(".", "_") + "\n")
                            else:
                                new_protein_handle.write(line)
            except IOError:
                raise Exception(f"Error reading or writing to file: {original_protein_file} or {protein_file}")
            # delete the original file
            os.remove(original_protein_file)
        # create checkpoint file
        try:
            with open(output.done, 'w') as done_handle:
                done_handle.write("done")
        except IOError:
            raise Exception(f"Error writing to file: {output.done}")


# to be honest, this should have been taken care of in the downloads commands
# but I only remembered to implement it later, and I do not want to rerun the entire download
# once more...
rule delete_ncbi_readme:
    input:
        genomic_download = "data/checkpoints_dataprep/{taxon}_A04_download.done",
        annotated_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl",
        blank_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl",
    output:
        done = "data/checkpoints_dataprep/{taxon}_A10_delete_ncbi_readme.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        # get all README files for a taxon
        df_anno = pd.read_csv(input.annotated_tbl_path, sep="\t", usecols=['species']) # here also proteins
        df_blank = pd.read_csv(input.blank_tbl_path, sep="\t", usecols=['species'])
        # concatenate the species of df_anno to df_blank
        df = pd.concat([df_anno, df_blank])
        genome_files = [f"data/species/{species.replace(' ', '_')}/README.md" for species in df_blank['species']]
        for genome_file in genome_files:
            # if file exists
            if os.path.exists(genome_file):
                os.remove(genome_file)
        # create checkpoint file
        try:
            with open(output.done, 'w') as done_handle:
                done_handle.write("done")
        except IOError:
            raise Exception(f"Error writing to file: {output.done}")


rule find_organelles:
    input:
        annotated_tbl_path = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl",
        download_done = "data/checkpoints_dataprep/{taxon}_A09_shorten_genomic_headers.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A09_b.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        logfile="data/checkpoints_dataprep/{params.taxon}_A09_b.log"
        echo "" > $logfile
        declare -a species_list
        while IFS=$'\t' read -r accession species others; do
            if [[ "$accession" == "accession" ]]; then
                continue
            fi
            modified_species="${{species// /_}}"
            species_list+=("$modified_species")
        done < {input.annotated_tbl_path}
        for species in "${{species_list[@]}}"; do
            bash scripts/find_organelles.sh data/species/${{species}}/genome/genome.fa > data/species/${{species}}/genome/organelles.lst
        done
        touch {output.done}
        """
        

# This rule performs a number of processing steps on the reference annotation file.
# The steps are necessary because we modified the FASTA headers, before.
# Among others, this rule extracts pseudogenes into a file pseudo.gff3, that is
# not directly needed as input for the DeepLearner, but for evaluation of 
# prediction accuracy, later.
rule select_pseudo:
    input:
        download_done = "data/checkpoints_dataprep/{taxon}_A09_b.done",
        annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A11_select_pseudo.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:devel"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"
        logfile="data/checkpoints_dataprep/{params.taxon}_A11_select_pseudo.log"
        echo "" > $logfile
        # Initialize the species_list as an array
        declare -a species_list
        # Read from the input file and modify the species names
        while IFS=$'\t' read -r accession species others; do
            if [[ "$accession" == "accession" ]]; then
                continue
            fi
            # Replace spaces with underscores in the species name
            modified_species="${{species// /_}}"
            species_list+=("$modified_species")
        done < {input.annotated_tbl}
        for species in "${{species_list[@]}}"; do
            echo "Processing species: ${{species}}" >> $logfile
            echo "grep '^>' data/species/${{species}}/genome/genome.fa > data/species/${{species}}/annot/deflines" >> $logfile
            grep '^>' data/species/${{species}}/genome/genome.fa > data/species/${{species}}/annot/deflines
            echo "grep -v -f data/species/${{species}}/genome/organelles.lst data/species/${{species}}/annot/deflines | cut -f1 -d' ' | cut -b2- > data/species/${{species}}/annot/z" >> $logfile
            grep -v -f data/species/${{species}}/genome/organelles.lst data/species/${{species}}/annot/deflines | cut -f1 -d' ' | cut -b2- > data/species/${{species}}/annot/z
            echo "paste data/species/${{species}}/annot/z data/species/${{species}}/annot/z > data/species/${{species}}/annot/list.tbl" >> $logfile
            paste data/species/${{species}}/annot/z data/species/${{species}}/annot/z > data/species/${{species}}/annot/list.tbl
            echo "rm data/species/${{species}}/annot/z data/species/${{species}}/annot/deflines" >> $logfile
            rm data/species/${{species}}/annot/z data/species/${{species}}/annot/deflines
            echo "sed ..." >> $logfile 
            sed -e 's/\ttmRNA\t/\tmRNA\t/g' -e 's/=tmRNA;/=mRNA;/g' data/species/${{species}}/annot/annot.gff3 | grep -v '\tinverted_repeat\t' > data/species/${{species}}/annot/annot_modified.gff3 && mv data/species/${{species}}/annot/annot_modified.gff3 data/species/${{species}}/annot/annot.gff3
            echo "gff_to_gff_subset.pl --in data/species/${{species}}/annot/annot.gff3 --out data/species/${{species}}/annot/tmp_annot.gff3 --list data/species/${{species}}/annot/list.tbl --col 2 --v &> /dev/null" >> $logfile
            gff_to_gff_subset.pl --in data/species/${{species}}/annot/annot.gff3 --out data/species/${{species}}/annot/tmp_annot.gff3 --list data/species/${{species}}/annot/list.tbl --col 2 --v &> /dev/null
            echo "##gff-version 3" > data/species/${{species}}/annot/annot.gff3
            echo "Running probuild..." >> $logfile
            probuild --stat_fasta --seq data/species/${{species}}/genome/genome.fa | cut -f1,2 | tr -d '>' | grep -v '^$' | awk '{{print "##sequence-region  " $1 "  1 " $2}}' >> data/species/$species/annot/annot.gff3
            cat data/species/${{species}}/annot/tmp_annot.gff3 | grep -v "^#" >> data/species/${{species}}/annot/annot.gff3
            echo "rm data/species/${{species}}/annot/tmp_annot.gff3" >> $logfile
            rm data/species/${{species}}/annot/tmp_annot.gff3
            echo "gt gff3 -force -tidy -sort -retainids -checkids -o data/species/${{species}}/annot/tmp_annot.gff3 data/species/${{species}}/annot/annot.gff3" >> $logfile
            gt gff3 -force -tidy -sort -retainids -checkids -o data/species/${{species}}/annot/tmp_annot.gff3 data/species/${{species}}/annot/annot.gff3
            echo "mv data/species/${{species}}/annot/tmp_annot.gff3 data/species/${{species}}/annot/annot.gff3" >> $logfile
            mv data/species/${{species}}/annot/tmp_annot.gff3 data/species/${{species}}/annot/annot.gff3
            echo "select_pseudo_from_nice_gff3.pl data/species/${{species}}/annot/annot.gff3 data/species/${{species}}/annot/pseudo.gff3" >> $logfile
            select_pseudo_from_nice_gff3.pl data/species/${{species}}/annot/annot.gff3 data/species/${{species}}/annot/pseudo.gff3
        done
        touch {output.done}
        """


# This rule uses AGAT to add introns. They are by default not contained in the
# reference annotation. AGAT also fixes some format inconsistency issues.
rule add_introns:
    input:
        download_done = "data/checkpoints_dataprep/{taxon}_A11_select_pseudo.done",
        annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A12_add_introns.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"
        logfile="data/checkpoints_dataprep/{params.taxon}_A12_add_introns.log"
        echo "" > $logfile
        # Initialize the species_list as an array
        declare -a species_list
        # Read from the input file and modify the species names
        while IFS=$'\t' read -r accession species others; do
            if [[ "$accession" == "accession" ]]; then
                continue
            fi
            # Replace spaces with underscores in the species name
            modified_species="${{species// /_}}"
            species_list+=("$modified_species")
        done < {input.annotated_tbl}
        for species in "${{species_list[@]}}"; do
            echo "Processing species: ${{species}}" >> $logfile
            ## ACHTUNG: FEHLERQUELLE MIT EMBL HIER:
            # Schritt 1 zur Loesung: grep -v -P "ID=id-[^;]+;Parent=gene-" annot.gff3 > tmp2_annot.gff3
            # ersetze annot.gff3 durch tmp2_annot.gff3
            grep -v -E "ID=id-[^;]+;Parent=gene-" data/species/${{species}}/annot/annot.gff3 > data/species/${{species}}/annot/tmp2_annot.gff3
            echo "agat_sp_add_introns.pl --gff data/species/${{species}}/annot/tmp2_annot.gff3 --out data/species/${{species}}/annot/annot_tmp.gff3" >> $logfile
            agat_sp_add_introns.pl --gff data/species/${{species}}/annot/tmp2_annot.gff3 --out data/species/${{species}}/annot/annot_tmp.gff3 &>> $logfile
        done
        touch {output.done}
        """

# This rule converts the modified gff3 to gtf format.
rule gff3_to_gtf:
    input:
        download_done = "data/checkpoints_dataprep/{taxon}_A12_add_introns.done",
        annotated_tbl = "data/checkpoints_dataprep/{taxon}_A03_annotated.tbl"
    output:
        done = "data/checkpoints_dataprep/{taxon}_A13_gtf.done"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:devel"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"
        logfile="data/checkpoints_dataprep/{params.taxon}_A13_gff3_to_gtf.log"
        echo "" > $logfile
        # Initialize the species_list as an array
        declare -a species_list
        # Read from the input file and modify the species names
        while IFS=$'\t' read -r accession species others; do
            if [[ "$accession" == "accession" ]]; then
                continue
            fi
            # Replace spaces with underscores in the species name
            modified_species="${{species// /_}}"
            species_list+=("$modified_species")
        done < {input.annotated_tbl}
        for species in "${{species_list[@]}}"; do
            echo "Processing species: ${{species}}" >> $logfile
            # remove lines containing tRNA to avoid the missing parent problem 
            grep -vi "trna" data/species/${{species}}/annot/annot_tmp.gff3 > data/species/${{species}}/annot/annot_tmp2.gff3 && mv data/species/${{species}}/annot/annot_tmp2.gff3 data/species/${{species}}/annot/annot_tmp.gff3
            echo "gff3_to_gtf.pl data/species/${{species}}/annot/annot_tmp.gff3 data/species/${{species}}/annot/annot.gtf" >> $logfile
            gff3_to_gtf.pl data/species/${{species}}/annot/annot_tmp.gff3 data/species/${{species}}/annot/annot.gtf &>> $logfile
            echo "rm data/species/${{species}}/annot/list.tbl data/species/${{species}}/annot/annot.gff3 data/species/${{species}}/annot/annot_tmp.gff3" >> $logfile
            rm data/species/${{species}}/annot/list.tbl data/species/${{species}}/annot/annot.gff3 data/species/${{species}}/annot/annot_tmp.gff3 data/species/${{species}}/annot/tmp2_annot.gff3
        done
        touch {output.done}
        """
