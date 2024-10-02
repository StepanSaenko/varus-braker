#!/bin/bash

OPTSTRING=":s:t:f:o:"
spid=""
threads=""
output_folder=""
output=""

while getopts ${OPTSTRING} opt; do
	case ${opt} in
		s) spid="${OPTARG}";;
		t) threads="${OPTARG}";;
		f) output_folder="${OPTARG}";;
		o) output="${OPTARG}";;
		?) echo "Invalid option"; exit 1;;
	esac
done 

wd=${PWD}
mkdir ${output_folder}
cp data/species/${spid}/genome/genome.fa ${output_folder}/genome.fa
cd ${output_folder}
log=${output_folder}/${spid}_masking.log
touch $log
echo "BuildDatabase -name ${spid} -dir data/species/${spid}/genome" &>> $log
BuildDatabase -name ${spid} genome.fa &>> $log
echo "RepeatModeler -database ${spid} -threads ${threads} -LTRStruct" &>> $log
RepeatModeler -database ${spid} -threads ${threads} -LTRStruct &>> $log
echo "RepeatMasker -pa ${threads} -xsmall -lib ${spid}-families.fa" &>> $log
RepeatMasker -pa ${threads} -xsmall -lib ${spid}-families.fa genome.fa &>> $log
cp genome.fa.masked $wd/data/species/${spid}/genome/genome.fa.masked
cp ${spid}-families.fa $wd/data/species/${spid}/genome/${spid}-families.fa
cp ${spid}-families.stk $wd/data/species/${spid}/genome/${spid}-families.stk
cp ${spid}_masking.log $wd/data/checkpoints_annotate/${spid}_masking.log
cd $wd
rm -rf ${output_folder}
touch ${output}
