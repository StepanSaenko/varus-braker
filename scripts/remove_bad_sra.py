#!/usr/bin/env python3

import re
import argparse
import os

parser = argparse.ArgumentParser(description='Remove bad SRA entries from a list of SRA accessions and cleanup the files')
parser.add_argument('-t', '--threshold', type=float, default=20, help='Alignment rate threshold (default: 20)')
parser.add_argument('-l', '--hisat2_log', type=str, required=True, help='Hisat2 log file')
parser.add_argument('-i', '--input', type=str, required=True, help='Input file with SRA accessions')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file with SRA accessions')
parser.add_argument('-p', '--path', type=str, required=True, help='Path to folders with fastq and hisat2 files for deletion')
args = parser.parse_args()

bad_sra = []
try:
    with open(args.hisat2_log, 'r') as f:
        for line in f:
            alignment_rate_float = None
            sra_id = re.search(r'\/([A-Z0-9]+)\.sam', line)
            if sra_id:
                sra_id_string = sra_id.group(1)
            # match the number before % 0.01% overall alignment rate                                                                                                                   
            alignment_rate = re.search(r'([^%]+)% overall alignment rate', line)
            if alignment_rate:
                alignment_rate_float = float(alignment_rate.group(1))
            # we must assume to have some metagenomic libraries, and they may have pretty bad but still useful reads,                                                                  
            # only remove the really bad ones, 20% is an arbitrary threshold for this                                                                                                  
            if alignment_rate_float is not None and alignment_rate_float < args.threshold:
                bad_sra.append(sra_id_string)
                print("I am appending SRA to bad list for " + str(sra_id_string))
except IOError:
            raise Exception(f"Error reading from: {args.hisat2_log}")

try:
    with open(args.input, 'r') as f:
        lines = f.readlines()
except IOError:
    raise Exception(f"Error reading from file: {args.input}")

no_rnaseq_anymore = []
try:
    with open(args.output, 'w') as f:
        for line in lines:
            line = line.strip()
            species, sra_ids = line.split('\t')
            # in species, replace the whitespaces by underscores
            species = species.replace(' ', '_')
            sra_ids = sra_ids.split(',')
            bad_sra_ids = sra_ids.copy()
            sra_ids = [sra_id for sra_id in sra_ids if sra_id not in bad_sra]
            bad_sra_ids = [sra_id for sra_id in bad_sra_ids if sra_id in bad_sra]
            if sra_ids:
                f.write(f"{species}\t{','.join(sra_ids)}\n")
            else:
                no_rnaseq_anymore.append(species)
except IOError:
    raise Exception(f"Error writing to file {args.output}!")

