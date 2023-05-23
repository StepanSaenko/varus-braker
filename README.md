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

See table_example.tbl
Columns are divided by tab
You may use without DNA-link, will try to find genome in ncbi-databse.
RNA-files should have "\_1" and "\_2" parts in names.

```
SPECIES_NAME DNA_LINK RNA_LINKS[optional,ID1_1,ID1_2,ID2_1,ID2_2 ... ]
Drosophila melanogaster https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz   /home/user/genomes/RNA/SRR13030903_1.fastq,/home/user/genomes/RNA/SRR13030903_2.fastq
Tribolium castaneum     https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.fna.gz 
Bombyx mori
#You may add commentaries into the table and also use local files, need absolute path
Drosophila Melanogaster /home/user/genomes/Drosophila/genome.fasta.masked
Populus trichocarpa     /home/user/genomes/Populus/genome.fasta.masked
Arabidopsis thaliana
Caenorhabditis elegans



varus-braker.py --input table.txt
```


## Configuration

File config.ini should have pathways to all executables. 
```
[VARUS]
varus_path = /home/saenkos/VARUS2/VARUS
hisat2_path = /home/saenkos/hisat2
sratoolkit_path = /home/saenkos/sratoolkit.3.0.2-ubuntu64/bin
batchsize = 100000
maxbatches = 1000

[BRAKER]
braker_cmd = /home/saenkos/braker3/BRAKER/scripts/braker.pl 
augustus_config_path = /home/saenkos/Augustus/config
augustus_bin_path = /home/saenkos/Augustus/bin
augustus_scripts_path = /home/saenkos/Augustus/scripts
diamond_path = /home/saenkos/
prothint_path = /home/saenkos/ProtHint/bin
genemark_path = /home/saenkos/GeneMark-ETP/bin
orthodb_path = /home/saenkos/orthodb-clades
excluded = species


[SLURM_ARGS]
partition = snowball
cpus_per_task = 36
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
Please choose the optimal batchsize and maxbatches values: huge batchsize could lead to increasing running time.

## Accuracies

|                                             | Gene  |       | Transcript |       | Exon  |       |
|:-------------------------------------------:|:-----:|:-----:|:----------:|:-----:|:-----:|:-----:|
|                                             |  Sn   |  Sp   |     Sn     |  Sp   |  Sn   |  Sp   |
|          *Caenorhabditis elegans*           |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 71.65 | 83.64 |   54.17    | 74.62 | 80.16 | 94.40 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 66.19 | 82.86 |   49.53    | 74.93 | 74.30 | 94.65 |
|  batchsize = 75000 maxbatches = 600 Total = 45M   | 71.54 | 85.34 |   54.21    | 76.23 | 78.67 | 94.84 |
|                                             |       |       |            |       |       |       |
|           *Arabidopsis thaliana*            |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 82.63 | 82.86 |   57.66    | 77.95 | 82.44 | 93.96 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 82.54 | 82.96 |   57.58    | 78.09 | 82.20 | 94.01 |
|  batchsize = 75000 maxbatches = 600 Total = 45M   | 82.64 | 83.08 |   57.64    | 78.14 | 82.12 | 94.08 |
|                                             |       |       |            |       |       |       |
|            *Populus trichocarpa*            |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 77.93 | 88.46 |   63.15    | 80.83 | 85.48 | 94.50 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 78.16 | 88.53 |   63.29    | 81.38 | 85.57 | 94.69 |
|  batchsize = 75000 maxbatches = 600 Total = 45M   | 78.33 | 88.62 |   63.40    | 81.52 | 85.59 | 94.82 |
|                                             |       |       |            |       |       |       |
|            *Medicago truncatula*            |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 50.28 | 72.97 |   50.28    | 65.27 | 77.79 | 88.07 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 50.53 | 72.91 |   50.53    | 65.26 | 77.79 | 88.22 |
|  batchsize = 75000 maxbatches = 600 Total = 45M   | 50.56 | 72.95 |   50.56    | 65.56 | 77.71 | 88.46 |
|                                             |       |       |            |       |       |       |
|         *Parasteatoda tepidariorum*         |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 49.02 | 59.55 |   41.56    | 53.45 | 55.76 | 86.67 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 46.47 | 59.71 |   39.07    | 53.89 | 51.97 | 87.24 |
|                                             |       |       |            |       |       |       |
|                *Danio rerio*                |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 56.67 | 69.55 |   35.42    | 66.57 | 66.07 | 91.98 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 56.79 | 70.72 |   35.36    | 67.76 | 65.34 | 92.60 |
|  batchsize = 75000 maxbatches = 600 Total = 45M   | 56.84 | 71.11 |   35.30    | 68.06 | 65.85 | 92.73 |
|                                             |       |       |            |       |       |       |
|          *Drosophila Melanogaster*          |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 86.45 | 90.11 |   60.03    | 82.61 | 83.88 | 95.34 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 86.89 | 89.85 |   60.07    | 81.99 | 83.47 | 94.95 |
|  batchsize = 75000 maxbatches = 600 Total = 45M   | 87.85 | 90.31 |   61.19    | 83.11 | 84.16 | 95.43 |
|                                             |       |       |            |       |       |       |
|           *Solanum lycopersicum*            |       |       |            |       |       |       |
| batchsize = 1500000 maxbatches = 200 Total = 300M | 46.34 | 47.79 |   37.54    | 47.50 | 83.88 | 94.52 |
| batchsize = 100000 maxbatches = 1000 Total = 100M | 46.23 | 47.76 |   37.43    | 47.42 | 83.40 | 94.47 |

