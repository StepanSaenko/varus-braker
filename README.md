# varus-braker
Small pipeline combining raw genome-data and VARUS.bam (or RNA-seq data) and run BRAKER3


# Software dependencies:


BRAKER3 from https://github.com/Gaius-Augustus/BRAKER

The best way to use singularity container:

```
singularity build braker3.sif docker://teambraker/braker3:latest
```

VARUS from https://github.com/Gaius-Augustus/VARUS

BUSCO from https://busco.ezlab.org/

hisat2 https://github.com/DaehwanKimLab/hisat2

sra-toolkit https://github.com/ncbi/sra-tools requeired version >= 3.0.2 
Important: GeneMark-ETP has sra-tools included, but v3.0.1

# Usage:

Input should be a table with at least 1 column:

```
See table_example.tbl
Columns are divided by space or tab

SPECIES_NAME DNA_LINK RNA_LINKS[optional,ID1_1,ID1_2,ID2_1,ID2_2 ... ]
Drosophila melanogaster https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz   /home/saenkos/homework/try/SRR13030903_1.fastq,/home/saenkos/homework/try/SRR13030903_2.fastq
Tribolium castaneum     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.fna.gz 
Bombyx mori


varus-braker.py --input table.txt [optional] --batchsize 100000 --maxBatches 5000 --clade arthropoda
```


# Configuration

File config.ini should have pathways to all executables. 
```
[VARUS]
varus_path = /home/saenkos/VARUS/VARUS2
hisat2_path = /home/saenkos/hisat2
sratoolkit_path = /home/saenkos/sratoolkit.3.0.2-ubuntu64/bin


[BRAKER]
braker_cmd = /home/saenkos/braker3/BRAKER/scripts/braker.pl 
augustus_config_path = /home/saenkos/Augustus/config
augustus_bin_path =
augustus_scripts_path = /home/saenkos/Augustus/scripts
diamond_path = /home/saenkos/
prothint_path = /home/saenkos/ProtHint/bin
genemark_path = /home/saenkos/GeneMark-ETP/bin

[SLURM_ARGS]
partition = snowball
cpus_per_task = 48
module_load = module load bedtools
```

In case using singularity containter :
```
[VARUS]
varus_path = /home/saenkos/VARUS/VARUS2
hisat2_path = /home/saenkos/hisat2
sratoolkit_path = /home/saenkos/sratoolkit.3.0.2-ubuntu64/bin


[BRAKER]
braker_cmd = singularity exec braker.sif braker.pl 
genemark_path = /home/saenkos/GeneMark-ETP/bin

[SLURM_ARGS]
partition = snowball
cpus_per_task = 48
module_load = module load bedtools
    module load singularity
