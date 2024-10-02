#!/usr/bin/env python3

import argparse
import re

parser = argparse.ArgumentParser(description='Build busco command')
parser.add_argument('-g', '--gff3', type=str, required=True, help='gff3 file')
args = parser.parse_args()

def check_for_gene_entry(lines, i, new_lines, seen_genes):
    #check for gene entry
    if re.match(r'.+\t.+\tgene', lines[i]):
        if lines[i] not in seen_genes:
            new_lines.append(lines[i])
            seen_genes.add(lines[i])
        return
    elif re.match(r'.+\t.+\t(CDS|exon|intron)', lines[i]) or re.match(r'.+\t.+\t.RNA', lines[i]):
        new_lines.append(lines[i])
        return
    #missing gene
    else:
        for j in range(i + 1, len(lines)):
            if re.match(r'.+\t.+\tgene', lines[j]):
                next_gene = j
                break
            else:
                next_gene = len(lines)
        for j in range(i + 1, next_gene):
            if re.match(r'.+\t.+\t.RNA', lines[j]):
                columns = lines[j].split('\t')
                name_match = re.search(r'ID=rna-([^;]+)', columns[8])
                if name_match:
                    name = name_match.group(1)
                    new_line = f'{columns[0]}\t{columns[1]}\tgene\t{columns[3]}\t{columns[4]}\t{columns[5]}\t{columns[6]}\t{columns[7]}\tID=gene-{name};Name={name};gbkey=Gene;gene_biotype=protein_coding;locus_tag={name}\n'
                    if new_line not in seen_genes:
                        new_lines.append(new_line)
                        seen_genes.add(new_line)
                return

def check_for_rna_entry(lines, i, new_lines, seen_rna):
    #check for rna entry
    if re.match(r'.+\t.+\t.RNA', lines[i]):
        if lines[i] not in seen_rna:
            new_lines.append(lines[i])
            seen_rna.add(lines[i])
        return
    elif re.match(r'.+\t.+\t(CDS|exon|intron|gene)', lines[i]):
        new_lines.append(lines[i])
        return
    #missing rna
    else:
        for j in range(i + 1, len(lines)):
            if re.match(r'.+\t.+\t.RNA', lines[j]):
                next_rna = j
                break
            else:
                next_rna = len(lines)
        for j in range(i + 1, next_rna):
            if re.match(r'.+\t.+\t(CDS|exon)', lines[j]):
                columns = lines[j].split('\t')
                name_match = re.search(r'ID=(cds|exon)-([^;]+)', columns[8])
                rna_type=re.search(r'gbkey=([^;]+)', columns[8])
                note = re.search(r'Note=([^;]+)', columns[8])
                if name_match:
                    name = name_match.group(2)
                    new_line = f'{columns[0]}\t{columns[1]}\t{rna_type.group(1)}\t{columns[3]}\t{columns[4]}\t{columns[5]}\t{columns[6]}\t{columns[7]}\tID=rna-{name};Parent=gene-{name};Note={note.group(1)};gbkey={rna_type.group(1)};locus_tag={name}\n'
                    if new_line not in seen_rna:
                        new_lines.append(new_line)
                        seen_rna.add(new_line)
                return

def main():
    new_lines = []
    seen_rna = set()
    seen_genes = set()
    current_name = ''
    count = 0

    # Read the GFF3 file
    try:
        with open(args.gff3, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                #first lines
                if re.match(r'^##', line):
                    new_lines.append(line)
                    continue
                #lines where rna or genes might be missing
                else: 
                    parts = line.split('\t')
                    name = re.search(r'ID=(.+)-([^;]+)', parts[8])
                    name_match = name.group(2)
                    # First entry
                    if count == 0:
                        current_name = name_match
                        count += 1
                        check_for_gene_entry(lines, i, new_lines, seen_genes)
                        check_for_rna_entry(lines, i, new_lines, seen_rna)
                    # Following entries
                    else:
                        if name_match != current_name:
                            check_for_gene_entry(lines, i, new_lines, seen_genes)
                            check_for_rna_entry(lines, i, new_lines, seen_rna)
                            current_name = name_match
                        else:
                            new_lines.append(line)

    except IOError:
        raise Exception(f"Error reading from: {args.gff3}")

    #write the updated lines back to the file
    try:
        with open(args.gff3, 'w') as f:
            f.writelines(new_lines)
    except IOError:
        raise Exception(f"Error writing to: {args.gff3}")

if __name__ == "__main__":
    main()
