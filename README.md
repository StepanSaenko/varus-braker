## varus-braker
Small pipeline combining raw genome-data and VARUS.bam (or RNA-seq data) and run BRAKER3


## Software dependencies


BRAKER3 from https://github.com/Gaius-Augustus/BRAKER

The best way to use singularity container:

```
singularity build braker3.sif docker://teambraker/braker3:latest
```

VARUS from https://github.com/Gaius-Augustus/VARUS

BUSCO from https://busco.ezlab.org/

hisat2 https://github.com/DaehwanKimLab/hisat2

sra-toolkit https://github.com/ncbi/sra-tools requeired version >= 3.0.2 \
**Important: GeneMark-ETP has sra-tools included, but v3.0.1**

GeneMark-ETP https://github.com/gatech-genemark/GeneMark-ETP

## Usage

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


## Configuration

File config.ini should have pathways to all executables. 
```
[VARUS]
varus_path = /home/saenkos/VARUS/VARUS2
hisat2_path = /home/saenkos/hisat2
sratoolkit_path = /home/saenkos/sratoolkit.3.0.2-ubuntu64/bin


[BRAKER]
braker_cmd = /home/saenkos/braker3/BRAKER/scripts/braker.pl 
augustus_config_path = /home/saenkos/Augustus/config
augustus_bin_path = /home/saenkos/Augustus/bin
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
```  

## Accuracies

|                                      | Gene  |       | Transcript |       | Exon  |       |
|:------------------------------------:|:-----:|:-----:|:----------:|:-----:|:-----:|:-----:|
|                                      |  Sn   |  Sp   |     Sn     |  Sp   |  Sn   |  Sp   |
|       *Caenorhabditis elegans*       |       |       |            |       |       |       |
|  batchsize = 5000 maxbatches = 100   | 30.66 | 84.12 |   22.54    | 78.66 | 26.61 | 94.92 |
|  batchsize = 50000 maxbatches = 100  | 51.08 | 84.74 |   38.33    | 76.42 | 53.01 | 94.99 |
|  batchsize = 50000 maxbatches = 500  | 62.11 | 84.09 |   46.72    | 74.74 | 69.05 | 94.67 |
| batchsize = 50000 maxbatches = 1000  | 63.75 | 83.15 |   47.94    | 73.76 | 71.83 | 94.33 |
| batchsize = 100000 maxbatches = 100  | 54.22 | 84.48 |   40.73    | 75.85 | 57.76 | 94.91 |
| batchsize = 100000 maxbatches = 500  | 62.89 | 83.25 |   47.33    | 74.10 | 70.83 | 94.50 |
| batchsize = 100000 maxbatches = 1000 | 65.00 | 82.49 |   48.90    | 73.12 | 73.58 | 94.07 |
| batchsize = 200000 maxbatches = 500  | 62.31 | 82.82 |   48.42    | 73.77 | 72.61 | 94.28 |
| batchsize = 200000 maxbatches = 1000 | 61.57 | 78.70 |   45.96    | 69.84 | 70.10 | 93.30 |
| batchsize = 150000 maxbatches = 2000 | 61.24 | 76.71 |   45.81    | 67.54 | 70.24 | 92.64 |
|                                      |       |       |            |       |       |       |
|        *Arabidopsis thaliana*        |       |       |            |       |       |       |
|  batchsize = 5000 maxbatches = 100   | 44.03 | 73.17 |   30.02    | 66.07 | 33.37 | 89.10 |
|  batchsize = 50000 maxbatches = 100  | 68.18 | 80.29 |   47.26    | 71.30 | 68.02 | 92.86 |
|  batchsize = 50000 maxbatches = 500  | 74.14 | 81.30 |   51.50    | 72.24 | 75.63 | 93.14 |
| batchsize = 50000 maxbatches = 1000  | 75.48 | 81.34 |   52.46    | 72.29 | 77.25 | 93.07 |
| batchsize = 100000 maxbatches = 100  | 69.19 | 80.57 |   47.93    | 71.72 | 69.90 | 93.12 |
| batchsize = 100000 maxbatches = 500  | 75.22 | 81.36 |   52.24    | 72.34 | 77.06 | 93.13 |
| batchsize = 100000 maxbatches = 1000 | 76.16 | 81.14 |   52.88    | 72.23 | 78.33 | 93.07 |
| batchsize = 200000 maxbatches = 500  | 75.81 | 81.01 |   52.66    | 72.28 | 77.96 | 93.00 |
| batchsize = 200000 maxbatches = 1000 | 76.41 | 80.76 |   53.08    | 72.01 | 78.68 | 92.84 |
| batchsize = 150000 maxbatches = 2000 | 77.03 | 80.79 |   53.48    | 71.84 | 79.41 | 92.74 |
