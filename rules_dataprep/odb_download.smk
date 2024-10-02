# Rule to download unique orthoDB partitions
# WARNING, DANGER! We hard coded here a concatenation of stramenopiles and viridiplantae partitions
# this makes sense for diatoms, but for other small clades, similar steps may be required!
# Rule to download unique orthoDB partitions and handle special case for stramenopiles
rule download_orthodb_partitions:
    output:
        fasta=f"{config['BRAKER']['orthodb_path']}/{{odb_partition}}.fa"
    params:
        url=lambda wildcards: f"https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/{wildcards.odb_partition}.fa.gz",
        orthodb_path=lambda wildcards: config['BRAKER']['orthodb_path'],
        special_url="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz"
    shell:
        """
        # check if directory exists
        if [ ! -d {params.orthodb_path} ]; then
            mkdir -p {params.orthodb_path}
        fi
        # Download the main partition
        curl -o {params.orthodb_path}/{wildcards.odb_partition}.fa.gz {params.url}
        # Check if the partition is stramenopiles and handle accordingly
        if [ "{wildcards.odb_partition}" = "Stramenopiles" ]; then
            # Download the Viridiplantae partition
            curl -o {params.orthodb_path}/Viridiplantae.fa.gz {params.special_url}
            # Concatenate Stramenopiles and Viridiplantae
            zcat {params.orthodb_path}/Stramenopiles.fa.gz {params.orthodb_path}/Viridiplantae.fa.gz > {params.orthodb_path}/Stramenopiles_combined.fa 2> log
            # Move the combined file to the original Stramenopiles output, replacing it
            mv {params.orthodb_path}/Stramenopiles_combined.fa {output.fasta}
            # Clean up the Viridiplantae file
            rm {params.orthodb_path}/Viridiplantae.fa.gz {params.orthodb_path}/Stramenopiles.fa.gz
        else
            # Move the single partition to the final output
            gunzip {params.orthodb_path}/{wildcards.odb_partition}.fa.gz
            mv {params.orthodb_path}/{wildcards.odb_partition}.fa {output.fasta}
        fi
        """


