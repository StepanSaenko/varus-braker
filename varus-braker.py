#!/usr/bin/env

import subprocess
import time
from time import sleep
import os
import urllib.request
import urllib.error
import gzip
import bz2
import tarfile
import zipfile
import rarfile
import re
import requests
import argparse
import shutil
import random
from multiprocessing import Pool, Semaphore
import configparser
from contextlib import closing


#defaults

parser = argparse.ArgumentParser(description='Run braker for the branch of species')
parser.add_argument('--input', type=str, help='path to the input file, table should have column with species names, column with links to DNA-data, [optional] column with links to RNA-data',  default="list.txt")
config = configparser.ConfigParser()
config.read('config.ini')
args = parser.parse_args()
input_file_path = args.input
partitition = config.get('SLURM_ARGS', 'partition')
clade_list = ['metazoa', 'vertebrata', 'viridiplantae', 'arthropoda', 'eukaryota', 'fungi', 'stramenopiles']
ortho_path = config.get('BRAKER', 'orthodb_path') + '/species/'
excluded = config.get('BRAKER', 'excluded')
JOB_LIST = []
main_dir = os.getcwd()

def remove_symbols_after_first_space(filename):
    lines = []
    with open(filename, "r") as file:
        for line in file:
            if line.startswith(">"):
                line = re.sub(" .*", "", line)
            lines.append(line)
    return lines

def decompress_file(file_path):
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rb") as gz_file:
            with open(file_path[:-3], "wb") as unzipped_file:
                unzipped_file.write(gz_file.read())
        file_path = file_path[:-3]
    elif file_path.endswith(".bz2") or file_path.endswith(".bzip2"):
        with bz2.open(file_path, "rb") as bz2_file:
            with open(file_path[:-4], "wb") as uncompressed_file:
                uncompressed_file.write(bz2_file.read())
        file_path = file_path[:-4]
    elif file_path.endswith(".tar"):
        with tarfile.open(file_path, "r") as tar:
            tar.extractall()
        file_path = file_path[:-4]
    elif file_path.endswith(".zip"):
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            zip_ref.extractall()
        file_path = file_path[:-4]
    elif file_path.endswith(".rar"):
        with rarfile.RarFile(file_path, "r") as rar_ref:
            rar_ref.extractall()
        file_path = file_path[:-4]
    else:
        print("Uncompressed file of wrong archive format: ", file_path)
        return (file_path)
    return file_path

def protein_data(species_name):
    """
    Choose the right protein data
    """
    subdirs = [f for f in os.listdir(ortho_path) if os.path.isdir(os.path.join(ortho_path, f))]
    if species_name in subdirs:
        protein_file = os.path.join(ortho_path, subdirs[subdirs.index(species_name)],excluded+"_excluded.fa")
        if os.path.exists(protein_file):
            return protein_file
    # Find protein database
    # Make a request to the NCBI ESearch API to get the taxon ID for the search term
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_params = {"db": "taxonomy", "term": species_name}
    delay = random.randint(1, 15)
    sleep(delay)
    esearch_response = requests.get(esearch_url, params=esearch_params)
    #print(esearch_response)
    if "429" in esearch_response:
        sleep(delay*2)
        esearch_response = requests.get(esearch_url, params=esearch_params)
    #print("2nd try:", esearch_response)     
    match_tax = re.search(r"<Id>(\d+)</Id>", esearch_response.text)
    if match_tax: 
        taxon_id = match_tax.group(1)
    else:
        taxon_id = '2759'
    # Make a request to the NCBI EFetch API to get the taxonomic lineage for the taxon ID
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    efetch_params = {"db": "taxonomy", "id": taxon_id, "retmode": "xml"}
    efetch_response = requests.get(efetch_url, params=efetch_params)
    lineage_init = re.search(r"<Lineage>(.*?)</Lineage>", efetch_response.text)
    if lineage_init:
        lineage = lineage_init.group(1)
    else: 
        lineage = 'cellular organisms; Eukaryota'
    #print("114: lineage :", lineage, "for ", species_name)

    clades = lineage.split("; ")
    for i in range(len(clades)-1, -1, -1):
        if clades[i].lower() in clade_list:
            #print("First matching clade from the right:", clades[i],  "for ", species_name)
            protein_file = clades[i].capitalize()+".fa"
            break
    print("protein file is: ", protein_file, "for ", species_name)
    if os.path.exists(protein_file):
        print(f"{protein_file} already exists in the directory.")
    else:
        subprocess.run(["wget", "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/"+clades[i].capitalize()+".fa.gz"])
        subprocess.run(["gzip", "-d", clades[i].capitalize()+".fa.gz"])
        protein_file = subprocess.run(["readlink", "-f", clades[i].capitalize()+".fa"], stdout=subprocess.PIPE).stdout.decode().strip()
    proteins_file_path = str(os.path.abspath(protein_file))
    # Return the path to the 'proteins.fasta' file
    return proteins_file_path

def rename_fasta(input_file):
    """
    Rename headers in a FASTA file.

    Parameters:
        input_file (str): The path to the input FASTA file.

    Returns:
        output_file (str): The path to the output FASTA file.
        translation_table (dict): A dictionary mapping old headers to new headers.
    """
    output_file = os.path.splitext(input_file)[0] + "_renamed" + os.path.splitext(input_file)[1]

    translation_table = {}

    with open(input_file, 'r') as in_f, open(output_file, 'w') as out_f:
        count = 0
        for line in in_f:
            if line.startswith('>'):
                count += 1
                old_header = line.strip()[1:]
                new_header = f'seq{count}'
                translation_table[old_header] = new_header
                out_f.write(f'>{new_header}\n')
            else:
                out_f.write(line)

    with open(input_file + '.translation_table.txt', 'w') as tt_f:
        for old_header, new_header in translation_table.items():
            tt_f.write(f'{old_header}\t{new_header}\n')

    return output_file, translation_table



omamerh5_file = "LUCA.h5"

if os.path.isfile(omamerh5_file):
    print(f"{omamerh5_file} already exists in the directory.")
else:
    subprocess.run(["wget", "https://omabrowser.org/All/LUCA.h5"])
    omamerh5_file = subprocess.run(["readlink", "-f", "LUCA.h5"], stdout=subprocess.PIPE).stdout.decode().strip()
omamerh5_file_path = str(os.path.abspath(omamerh5_file))
# Print the path to the 'proteins.fasta' file
#print("Path to LUCA.h5: ", omamerh5_file_path)

def up_low(path):

    # Define a regular expression pattern to match upper-case DNA strings
    UPPER_PATTERN = r'^[ACGTN\s]+$'
    # Define a regular expression pattern to match lower-case DNA strings
    LOWER_PATTERN = r'^[acgtn\s]+$'
    # Open the fasta file
    #print("uplow path", path)
    # Read the contents of the file
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            fasta_contents = f.read() 
        #fasta_contents = f.read()
            if re.match(UPPER_PATTERN, fasta_contents):
                return("upper")
            elif re.match(LOWER_PATTERN, fasta_contents):
                return("lower")
            else:
                return("mixed")
    except UnicodeDecodeError as e:
        print("Error reading the file:", str(e))
        return("error")

def repeatmasking(dna_path, genus):

    repmask_script = """#!/bin/bash
#SBATCH --output=slurm-%j.out
#SBATCH -J RM
#SBATCH --get-user-env
#SBATCH -N 1 # number of nodes
#SBATCH -n 50
#SBATCH -p {}

###-LOAD PRE-INSTALLED MODULES

###-UPDATE YOUR INPUTS
INPUT={}
DB={}

###-SETUP DETECTION OF REPEATS
echo "Build database ..."
time BuildDatabase -name $DB $INPUT

echo "Run RepeatModeler ..."
time RepeatModeler -database $DB -threads 50

ln -s RM_*/consensi.fa.classified ./

echo "Run RepeatMasker ..."
time RepeatMasker -pa 50 -gff -lib consensi.fa.classified $INPUT
""".format(partitition, dna_path, genus)
    #print("RM script = :", repmask_script)
    repmask_process = subprocess.Popen(['sbatch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    repmask_process.stdin.write(repmask_script.encode())
    output_r, error_r = repmask_process.communicate()
    repmask_job_id = re.search(r"\d+", output_r.decode().strip()).group()
    #print("RepeatMasker Job ID:", repmask_job_id)
    dna_masked_path = dna_path + ".masked"
    #print("masked path : ", dna_masked_path)
    repmask_process.stdin.close()
    return(repmask_job_id, dna_masked_path)

def varus_run(dna_path, genus, species):
    current_dir = os.getcwd()
    #print("VARUS Current directory:", current_dir)
    new_dir = os.path.dirname(dna_path)
    os.chdir(new_dir)
    # get new current working directory
    new_current_dir = os.getcwd()
    name_id = str(genus) + '_' + str(species)
    if os.path.exists(f"{new_current_dir}/varus.err"):
        # Delete the file
        os.remove(f"{new_current_dir}/varus.err")
        #print("File varus.err has been deleted")
    #else:
        #print("File varus.err does not exist")
    if os.path.exists(f"{new_current_dir}/{name_id}/VARUS.bam"):
        os.chdir(current_dir)
        print("varus.bam is already exists")
        return ("1", f"{new_current_dir}/{name_id}/VARUS.bam")

    #parse config file
    varus_path = config.get('VARUS', 'varus_path')
    hisat2_path = config.get('VARUS', 'hisat2_path')
    sratoolkit_path = config.get('VARUS', 'sratoolkit_path')
    batchsize = str(config.get('VARUS', 'batchsize'))
    maxbatches = str(config.get('VARUS', 'maxbatches'))
    
    hisat2_export_path = f"export PATH={hisat2_path}:$PATH" if hisat2_path else ""
    sratoolkit_export_path = f"export PATH={sratoolkit_path}:$PATH" if sratoolkit_path else ""

    with open('VARUSparameters.txt', 'w') as f:
        f.write(f"--batchSize {batchsize}\n")
        f.write(f"--blockSize 5000\n")
        f.write(f"--components 1\n")
        f.write(f"--cost 0.001\n")
        f.write(f"--deleteLater 0\n")
        f.write(f"--estimator 2\n")
        f.write(f"--exportObservationsToFile 1\n")
        f.write(f"--exportParametersToFile 1\n")
        f.write(f"--fastqDumpCall fastq-dump\n")
        f.write(f"--genomeDir ./genome/\n")
        f.write(f"--lambda 10.0\n")
        f.write(f"--lessInfo 1\n")
        f.write(f"--loadAllOnce 0\n")
        f.write(f"--maxBatches {maxbatches}\n")
        f.write(f"--mergeThreshold 10\n")
        f.write(f"--outFileNamePrefix ./\n")
        f.write(f"--pathToParameters ./VARUSparameters.txt\n")
        f.write(f"--pathToRuns ./\n")
        f.write(f"--pathToVARUS {varus_path}/Implementation\n")
        f.write(f"--profitCondition 0\n")
        f.write(f"--pseudoCount 1\n")
        f.write(f"--qualityThreshold 5\n")
        f.write(f"--randomSeed 1\n")
        f.write(f"--readParametersFromFile 1\n")
        f.write(f"--runThreadN 32\n")
        f.write(f"--verbosityDebug 1\n")

      
    slurm_varus = f"""#!/bin/bash
#SBATCH -o varus.%j.%N.out
#SBATCH -e varus.%j.%N.err
#SBATCH -J varus
#SBATCH --get-user-env
#SBATCH -N 1 # number of nodes
#SBATCH -n 36
#SBATCH -p {partitition}
{hisat2_export_path}
{sratoolkit_export_path}
{varus_path}/runVARUS.pl --aligner=HISAT --readFromTable=0 --createindex=1 --latinGenus={genus} --latinSpecies={species} --speciesGenome={os.path.basename(dna_path)} --logfile=varus.log 2>varus.err"""
    while True:
        process = subprocess.Popen(['sbatch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.stdin.write(slurm_varus.encode())
        output, error = process.communicate()
        job_id = re.search(r"\d+", output.decode().strip()).group()
        #print("VARUS Job ID:", job_id)
        process.stdin.close()
        varus_bam = f"{new_current_dir}/{name_id}/VARUS.bam"
        varus_err = f"{new_current_dir}/varus.err"
        #print("Varus.bam = ", varus_bam)
        #print("Varus.err = ", varus_err)
        time.sleep(10)
        while not os.path.exists(varus_err):
            time.sleep(30)
            continue
        time.sleep(240)
        with open(varus_err, 'r') as f:
            content = f.read()          
        if "error 200" in content or "error 500" in content or "error 429" in content:
                #if len(lines) >= 5 and ("500 Internal Server Error" in lines[3] or "ERROR 500: Internal Server Error." in lines[4] or "ERROR 429" in lines[4]):
            print("Server error 500/429. Waiting 300 seconds and trying again...")
            subprocess.run(['scancel ', job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            time.sleep(300)
            continue
        break
    os.chdir(current_dir)
    return(job_id, varus_bam)


def braker_run(dna_path, rna_path, genus, species, proteins_file_path):
    #print("braker_parameters : 1: ",dna_path,"2: ",genus,species,"3: ",rna_path, "->",os.path.basename(rna_path))
    current_dir = os.getcwd()
    new_dir = os.path.dirname(os.path.abspath(dna_path))
    os.chdir(new_dir)
    w_dir = genus+"_"+species+"_braker"
    # get new current working directory
    
    if os.path.exists(f"{os.getcwd()}/{w_dir}/braker.gtf"):
        os.chdir(current_dir)
        print("BRAKER.gff is already exists")
        return ("1", f"{os.getcwd()}/{w_dir}/braker.gtf")
    
    if os.path.basename(rna_path) == "NNNN":
        #run BRAKER2
        rna_subline = ""
    elif os.path.basename(rna_path)[-3:] == "bam":
        rna_subline = " --bam="+rna_path
    else:
        rna_paths = rna_path.split(',')
        prefixes = list(set([os.path.basename(s).split('_')[0] for s in rna_paths])) 
        rna_names  = ",".join(prefixes)
        rna_subline =  " --rnaseq_sets_ids="+rna_names+ " --rnaseq_sets_dirs="+os.path.dirname(os.path.abspath(rna_paths[0]))

    augustus_bin_path = config.get('BRAKER', 'augustus_bin_path')
    augustus_config_path = config.get('BRAKER', 'augustus_config_path')
    augustus_scripts_path = config.get('BRAKER', 'augustus_scripts_path')
    module_load = config.get('SLURM_ARGS', 'module_load')
    braker_cmd = config.get('BRAKER', 'braker_cmd')
    diamond_path = config.get('BRAKER', 'diamond_path')
    prothint_path = config.get('BRAKER', 'prothint_path')
    genemark_path = config.get('BRAKER', 'genemark_path')
    
    augustus_bin_arg = f"--AUGUSTUS_BIN_PATH={augustus_bin_path}" if augustus_bin_path else ''
    augustus_config_arg = f"--AUGUSTUS_CONFIG_PATH={augustus_config_path}" if augustus_config_path else ''
    augustus_scripts_arg = f"--AUGUSTUS_SCRIPTS_PATH={augustus_scripts_path}" if augustus_scripts_path else ''
    diamond_arg = f"--DIAMOND_PATH={diamond_path}" if diamond_path else ''
    prothint_arg = f"--PROTHINT_PATH={prothint_path}" if prothint_path else ''
    genemark_arg = f"--GENEMARK_PATH={genemark_path}" if genemark_path else ''
    genemark_export = f"export PATH={genemark_path}/tools:$PATH" if genemark_path else ''
    module_arg = f"{module_load}" if module_load else ''
    id = random.randint(100, 10000)
    slurm_braker = f"""#!/bin/bash
#SBATCH -o braker.%j.%N.out
#SBATCH -e braker.%j.%N.err
#SBATCH -J braker3
#SBATCH --get-user-env
#SBATCH -N 1 # number of nodes
#SBATCH -n 36
#SBATCH -p {partitition}

{genemark_export}
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


mkdir /tmp/saenkos-{id}
cp -r ./ /tmp/saenkos-{id}/
cp {proteins_file_path} /tmp/saenkos-{id}/
cd /tmp/saenkos-{id}/

{module_arg}
{braker_cmd} {augustus_config_arg} {augustus_bin_arg} {augustus_scripts_arg} \
{diamond_arg} {prothint_arg} --softmasking --useexisting {genemark_arg} --threads 12 \
--species={genus}_{species} --workingdir=./{w_dir} --prot_seq={proteins_file_path} --genome=./{os.path.basename(dna_path)} {rna_subline}

cd -
mv /tmp/saenkos-{id}/{w_dir} ./
rm -rf /tmp/saenkos-{id}
"""
    #print(slurm_braker)
    process = subprocess.Popen(['sbatch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.stdin.write(slurm_braker.encode())
    output, error = process.communicate()
    job_id = re.search(r"\d+", output.decode().strip()).group()
    #print("BRAKER3 Job ID:", job_id)
    process.stdin.close()
    gtf_file = f"{os.getcwd()}/{w_dir}/braker.gtf"
    print("braker.gtf is: ", gtf_file)
    os.chdir(current_dir)
    return(job_id,gtf_file)

def busco_run(braker_dir, out_name):
    slurm_busco = f"""#!/bin/bash
#SBATCH -o busco.%j.%N.out
#SBATCH -e busco.%j.%N.err
#SBATCH -J busco
#SBATCH --get-user-env
#SBATCH -N 1 # number of nodes
#SBATCH -n 36
#SBATCH -p {partitition}
    
    
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export NUMEXPR_MAX_THREADS=48

cd {braker_dir}
busco  -i braker.codingseq -c 16 -m geno -f --out ./BUSCO/{out_name} --auto-lineage-euk

"""
    process = subprocess.Popen(['sbatch'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.stdin.write(slurm_busco.encode())
    output, error = process.communicate()
    busco_job_id = re.search(r"\d+", output.decode().strip()).group()
    
    process.stdin.close()
    return(busco_job_id)

def process_line(line):
    RESULTS_FILE = "results.log"
    ERROR_FILE = "errors.log"
    rna_paths = []
    RM_failed = False
    dna_fail = False
    rna_fail = False
    dna_count = 0
    line_parts = re.split(r'\s+|\t+',str(line))
    genus = line_parts[0]
    species = line_parts[1]
    links = line_parts[2:]
    current_dir = os.getcwd()
    name_id = str(genus) + '_' + str(species)
    #print("Current directory:", current_dir, " for species: ", name_id )
    JOB_LIST.append(name_id) 

    protein_file_path = protein_data(name_id)
    
    # create a directory with the name
    os.makedirs(name_id, exist_ok=True)
    results_path = os.path.join(name_id, RESULTS_FILE)
    try:
        with open(results_path, 'a') as f:
            pass # File exists, do nothing
    except FileNotFoundError:
        with open(results_path, 'w') as f:
            pass # File doesn't exist, create it
    error_path = os.path.join(name_id, ERROR_FILE)
    try:
        with open(error_path, 'a') as f:
            pass # File exists, do nothing
    except FileNotFoundError:
        with open(error_path, 'w') as f:
            pass
    #print("links to files :",links)
    if len(links) > 0:
        if os.path.isabs(links[0]):
            print("It is a local DNA-file!")
            destination_file = os.path.join(os.getcwd(), name_id, os.path.basename(links[0]))
            if not os.path.exists(destination_file):
                shutil.copy(links[0], os.path.join(os.getcwd(), name_id))
            shutil.copy(links[0], destination_file)
            # Create a new pathway to the copied file in the directory
            dna_path = destination_file
            # Print the new file pathway
            #print(dna_path)
        else:
            print("downloading DNA data :")
            #filename = os.path.basename(links[0])
            try:
                response = requests.get(links[0], allow_redirects=True, timeout=300)
                               
                # Retrieve the filename from the response headers
                filename = os.path.basename(response.url)
                response.raise_for_status()
                dna_path = os.path.join(name_id, filename.split("?")[0])
                if not os.path.exists(dna_path):
                    with open(dna_path, "wb") as file:
                        file.write(response.content)
                #urllib.request.urlretrieve(links[0], f"{name_id}/{os.path.basename(links[0])}")
            except requests.Timeout:
                print("Request timed out. Aborting.")
                dna_fail = True
            except requests.HTTPError as e:
                print("HTTP error occurred:", str(e))
                dna_fail = True
            except requests.RequestException as e:
                print("An error occurred:", str(e))
                dna_fail = True
            else:
                dna_path = f"{name_id}/{os.path.basename(links[0])}"
        if not dna_fail:
            dna_path = decompress_file(dna_path)
            
    if len(links) == 0 or dna_fail:
        print(name_id, " is empty or wrong link, trying to use ncbi-datasets")
        os.chdir(name_id)
        organism = str(genus) + ' ' + str(species)
        while True:
            try:
                completed_process = subprocess.run(
                    ["datasets", "download", "genome", "taxon", organism, "--reference"],
                    capture_output=True,
                    text=True,  # This is important to get the output as text
                    check=True   # This will raise a CalledProcessError if the command returns a non-zero exit status
                )
            except subprocess.CalledProcessError as e:
                error_message = e.stderr.strip()  # Capturing the error message
                if "Error: No assemblies found that match selection" in error_message:
                    print("No matching assemblies found.")
                    os.chdir(current_dir)
                    JOB_LIST.remove(name_id) 
                    return("NCBI datasets fail.")
                elif "Error: Internal error (invalid zip archive). Please try again" in error_message:
                    print("Internal error with zip archive. Retrying...")
                    continue  # Retry the subprocess
                else:
                    print("An error occurred:", error_message)
                    os.chdir(current_dir)
                    JOB_LIST.remove(name_id) 
                    return("NCBI datasets fail.")
            break
        time.sleep(20)
        zip_filename = "ncbi_dataset.zip"
        with zipfile.ZipFile(zip_filename) as zip_file:
        # Loop through each file in the archive
            for filename in zip_file.namelist():
            # If the file has a ".fna" extension, move it to the current directory and break out of the loop
                if os.path.splitext(filename)[1] in [".fna", ".fasta", ".fa"]:
                # Extract the file to a temporary directory
                    zip_file.extract(filename, path = "temp")
                    destination_file_path = os.path.join(os.getcwd(), os.path.basename(filename))
                    # Move the file from the temporary directory to the current directory
                    if not os.path.exists(destination_file_path):
                        shutil.move(os.path.join("temp", filename), ".")
                        moved_file_path = os.path.join(os.path.basename(os.getcwd()), os.path.basename(filename))
                        #print(f"File found and moved: {moved_file_path}")
                        dna_count = dna_count + 1
                        break
                    else:
                        moved_file_path = destination_file_path
                        print(f"File {filename} already exists in the current directory.")
                        dna_count = dna_count + 1
                        break
        """
        with zipfile.ZipFile(zip_filename) as zip_file:
        # Loop through each file in the archive
            for filename in zip_file.namelist():
            # If the file has a ".fna" extension, move it to the current directory and break out of the loop
                if os.path.splitext(filename)[1] == ".fna" or os.path.splitext(filename)[1] == ".fasta" or os.path.splitext(filename)[1] == ".fa":
                # Extract the file to a temporary directory
                    zip_file.extract(filename, path="temp")
                    # Move the file from the temporary directory to the current directory
                    if not os.path.exists(os.path.join("temp", filename)):
                        shutil.move(os.path.join("temp", filename), ".")
                    moved_file_path = os.path.join(os.path.basename(os.getcwd()), os.path.basename(filename))
                    #print(f"File found and moved: {moved_file_path}")
                    dna_count = dna_count + 1
                    break"""
        shutil.rmtree("temp")
        dna_path = moved_file_path
        os.chdir(current_dir)                
       
        if dna_count == 0:
            print ("No DNA in link and in ncbi, check the link or  genus and species.")
            with open(error_path, 'a') as f:
                f.write(f"{name_id} : lost DNA data\n")
            os.chdir(current_dir)
            JOB_LIST.remove(name_id) 
            return("DNA fail.")          
        #RNA processing        
    with open(results_path, 'a') as f:
        f.write(f"Species: {name_id} : {dna_path}\n")
                
    if len(links) > 1:
        rna_links = links[1].split(',')
        rna_paths = []
        for link in rna_links:
            #check if it is a local file or file on a web server
            if os.path.isabs(link):
                rna_path = link
                rna_paths.append(rna_path)
                print("It is a local RNA-file!")
            else:
                print("downloading RNA data :")
                try:
                    urllib.request.urlretrieve(link, f"{name_id}/{os.path.basename(link)}")
                except urllib.error.HTTPError as e:
                    if e.code == 404:
                        print("RNA file not found, check the link. Will use VARUS instead.")
                        with open(error_path, 'a') as f:
                            f.write(f"{name_id} : lost RNA data\n")
                else:
                    rna_path = f"{name_id}/{os.path.basename(link)}"
                    rna_path = decompress_file(rna_path)
                    rna_paths.append(rna_path)
            #print("rna paths :", rna_paths)
            with open(results_path, 'a') as f:
                f.write(f"{rna_paths}\n")
    if up_low(dna_path) == "error":
        print("Something wrong with fasta file : ", name_id, " Check the error file...")
        with open(error_path, 'a') as f:
            f.write(f"{name_id} : something wrong with fasta\n")
        return 0
    elif up_low(dna_path) == "mixed":
        print("Mixed case in DNA file. No need to run RepeatMasker...")
        masked_dna_path = dna_path
    else:
        print("running RepeatMasker for ... ", name_id)
        RM_job_id, masked_dna_path = repeatmasking(dna_path, genus)
        while True:
            result = subprocess.run(['squeue', '-j', RM_job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            time.sleep(5000)
            output = result.stdout + result.stderr
            if RM_job_id not in output:
                if os.path.isabs(masked_dna_path):
                    if os.path.exists(masked_dna_path):
                        print("RepeatMasking is done. Masked file exists")
                        with open(results_path, 'a') as f:
                            f.write(f"{name_id} : RepeatMasking finished successfully wtih {RM_job_id}. {masked_dna_path}\n")
                        break
                    else:
                        RM_failed = True
                        with open(error_path, 'a') as f:
                            f.write(f"{name_id} : RepeatMasking fail with {RM_job_id}\n")
                        return("RepeatMasking fail")

    short_header_dna, tr_table = rename_fasta(masked_dna_path)
    if rna_paths:
        print("RNA files are presented. No need to run VARUS...")
        rna_file = ','.join(map(str, rna_paths))
    else:
        varus_job_id, rna_file = varus_run(short_header_dna, genus, species)
        print("VARUS is runnig with ", varus_job_id, ", .bam file will be: ", rna_file)
        time.sleep(10)
        while True:
            var_result = subprocess.run(['squeue', '-j', varus_job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            var_output = var_result.stdout + var_result.stderr
            if varus_job_id not in var_output:
                if os.path.isabs(rna_file):
                    if os.path.exists(rna_file):
                        print("VARUS is done. BAM file exists")
                        with open(results_path, 'a') as f:
                            f.write(f"{name_id} : VARUS finished successfully with {varus_job_id}. {rna_file}\n")
                        break
                    else:
                        print("VARUS has failed. BAM file does not exist. Check the error file...")
                        with open(error_path, 'a') as f:
                            f.write(f"{name_id} : VARUS failed with {varus_job_id}\n")
                        varus_failed = True
                        rna_file = "NNNN"
                        print("Will run BRAKER without RNA data. BRAKER2")
                        break
                else:
                    print("VARUS has failed. BAM file does not exist. Check the error file...")
                    with open(error_path, 'a') as f:
                        f.write(f"{name_id} : VARUS failed with {varus_job_id}\n")
                    varus_failed = True
                    rna_file = "NNNN"
                    print("Will run BRAKER without RNA data. BRAKER2")
                    break
            time.sleep(600)            
    print("rna_file :", rna_file)
    #print("____________________________________________________")
    print("BRAKER run for ", genus, species)
    braker_job_id, gtf_file = braker_run(short_header_dna, rna_file, genus, species, protein_file_path)
    print("BRAKER3 job =", braker_job_id)
    
    while True:
        br_result = subprocess.run(['squeue', '-j', braker_job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
       
        br_output = br_result.stdout + br_result.stderr
        if braker_job_id not in br_output:
            if os.path.isabs(gtf_file):
                if os.path.exists(gtf_file):
                    print("BRAKER is done. GTF file exists")
                    with open(results_path, 'a') as f:
                        f.write(f"{name_id} : BRAKER finished successfully with {braker_job_id}. {gtf_file}\n")
                        #f.write("____________________________________________________\n")
                    print("check the quality using files braker.aa (protein) and braker.codingseq(codingseq)")
                    #omark_dir = name_id+'_OMArk'
                    #omark_command = f"omark.py -f '{name_id}/{name_id}_braker/braker.codingseq' -d LUCA.h5 -o '{name_id}/{omark_dir}'"
                    #subprocess.run(omark_command, stdout=subprocess.PIPE, shell=True)
                    """We need to wait ~ a week until OMArk will update for python 3.11"""
                    bsc_job_id = busco_run(f"./{name_id}/{name_id}_braker/", name_id)
                    with open(results_path, 'a') as f:
                        f.write(f"{name_id} : BUSCO lauhcned with {bsc_job_id}\n")
                    #busco_dir = name_id+'_BUSCO'
                    #busco_command = f"busco  -i ./{name_id}/{name_id}_braker/braker.codingseq -c 16 -m geno -f --out ./{name_id}/{busco_dir} --auto-lineage-euk"
                    #print(busco_command)
                    #subprocess.run(busco_command, stdout=subprocess.PIPE, shell=True)
                    break
                else:
                    print("BRAKER has failed. GTF file does not exist. Check the error file...")
                    with open(error_path, 'a') as f:
                        f.write(f"{name_id} : BRAKER failed with {braker_job_id}\n")
                        f.write("________________________________________________\n")
                    braker_failed = True
                    JOB_LIST.remove(name_id) 
                    return("BRAKER fail")
                    break
        time.sleep(6000)

def process_with_semaphore(part):
    with semaphore:
        sleep(35)
        process_line(part)
        
if __name__ == '__main__':
    MAX_PARALLEL_JOBS = 25
    semaphore = Semaphore(MAX_PARALLEL_JOBS)
    with open(input_file_path, 'r') as f:
        next(f)
        lines = [line.strip() for line in f.readlines()]
    #print(lines)
    processed_lines = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            processed_lines.append(line)    
    #print(processed_lines)
    
    parts = [processed_lines]
    #parts = [lines[i:i+10] for i in range(0, len(lines), 10)]
    #print (parts)            
    with Pool() as pool:
        for part in processed_lines:
            sleep(35)
            pool.apply_async(process_with_semaphore, (part,))
            #pool.map(process_line, part)
    # Close the pool and wait for all processes to finish
        pool.close()
        pool.join()

    # Save the results to an output file

    try:
        with open('results_sum.log', 'a') as f:
            pass # File exists, do nothing
    except FileNotFoundError:
        with open('results_sum.log', 'w') as f:
            pass # File doesn't exist, create it
    try:
        with open('errors_sum.log', 'a') as f:
            pass # File exists, do nothing
    except FileNotFoundError:
        with open('errors_sum.log', 'w') as f:
            pass
    print("combinig results and errors files")
    with open("results_sum.log", "wb") as outfile, open("errors_sum.log", "wb") as errfile:
        # Walk through all directories in the main directory
        for dir_name in os.listdir(main_dir):
            dir_path = os.path.join(main_dir, dir_name)
            if os.path.isdir(dir_path):
            # Check if the current directory has the necessary files
                if "results.log" in os.listdir(dir_path):
                    # If it does, add the output to the OUTPUT.out file
                    with open(os.path.join(dir_path, "results.log"), "rb") as infile:
                        shutil.copyfileobj(infile, outfile)
                if "errors.log" in os.listdir(dir_path):
                    # If it does, add the output to the ERROR.out file
                    with open(os.path.join(dir_path, "errors.log"), "rb") as infile:
                        shutil.copyfileobj(infile, errfile)

    print("All done! Please check the results_sum.log and errors_sum.log files")
