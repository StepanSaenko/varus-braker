#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Build busco command')
parser.add_argument('-c', '--csv', type=str, required=True, help='csv file with all info about all species')
parser.add_argument('-s', '--species', type=str, required=True, help='Species name in data structure')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file name with annotation command')
#parser.add_argument('-b', '--braker', type=str, required=True, help='Braker file to be analyze')
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use')
parser.add_argument('-w', '--workingdir', type=str, default='.', help='Working directory')
parser.add_argument('-l', '--logfile', type=str, default='braker.log', help='Log file name')
args = parser.parse_args()

# read csv file
df = pd.read_csv(args.csv)

# find the row where column "species" has the value of args.species
row = df.loc[df['species'] == args.species]

cmd_g = ""
cmd_g += 'busco --in=data/species/' + str(args.species) + '/genome/genome.fa'  
cmd_g += ' --mode=genome'
cmd_g += ' --lineage_dataset=' + row['busco_lineage'].values[0]
cmd_g += ' --cpu=' + str(args.threads)
cmd_g += ' --out=/data/species/' + str(args.species) + '/busco/genomefile'
cmd_g += ' &> ' + args.logfile

cmd_b = ""
cmd_a = ""

if pd.isna(row['annotation_file'].values[0]):
    cmd_b += 'busco --in=data/species/' + str(args.species) + '/braker/braker.aa'
    cmd_b += ' --mode=proteins'
    cmd_b += ' --lineage_dataset=' + row['busco_lineage'].values[0]
    cmd_b += ' --cpu=' + str(args.threads)
    cmd_b += ' --out=/data/species/' + str(args.species) + '/busco/brakerfile'
    cmd_b += ' &> ' + args.logfile
else:
    cmd_a += 'busco --in=data/species/' + str(args.species) + '/prot/protein.faa'
    cmd_a += ' --mode=proteins'
    cmd_a += ' --lineage_dataset=' + row['busco_lineage'].values[0]
    cmd_a += ' --cpu=' + str(args.threads)
    cmd_a += ' --out=/data/species/' + str(args.species) + '/busco/annotfile'
    cmd_a += ' &> ' + args.logfile

# write command to output file
try:
    with open(args.output, 'w') as f:
        f.write(cmd_g)
        f.write('\n')
        if(len(cmd_a)>0):
            f.write(cmd_a)
            f.write('\n')
        if(len(cmd_b)>0):
            f.write(cmd_b)
            f.write('\n')

except IOError:
    print("Error: can't write to file" + args.output)
    exit(1)
