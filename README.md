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

# Usage:

Input should be a table with at least 1 column:

```
SPECIES_NAME DNA_LINK RNA_LINK(s)[optional]

```
# Configuration

File config.ini should have pathways to all executables. 
```
[VARUS]
varus_path = /home/saenkos/VARUS/VARUS2
hisat2_path = 
sratoolkit_path =

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
hisat2_path = 
sratoolkit_path =

[BRAKER]
braker_cmd = singularity exec braker.sif braker.pl 
genemark_path = /home/saenkos/GeneMark-ETP/bin

[SLURM_ARGS]
partition = snowball
cpus_per_task = 48
module_load = module load bedtools
    module load singularity
