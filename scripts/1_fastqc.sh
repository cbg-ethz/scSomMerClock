#!/bin/sh

module purge
module load gcccore/6.4.0 cutadapt/1.18-python-3.7.0

##$1: Sample name
##$2: Reference genome file
##$3: String identifyping the WGA protocol
sample=$1
REF=$2
WGA_LIBRARY=$3

cutadapt \
	--minimum-length 70 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    -o Processing/${sample}.trimmed_1.fastq.gz \
    -p Processing/${sample}.trimmed_2.fastq.gz \
    Raw_Data/${sample}_1.fastq.gz Raw_Data/${sample}_2.fastq.gz \
    > ./slurm_files/${sample}_Cutadapt.txt


if [[ ${WGA_LIBRARY} == "AMPLI-1" ]]
then
        Adapter1="GCTGTCAGTTAA"
        Adapter2="TTAACTGACAGCAGGAATCCCACT"
        adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2}"
elif [[ ${WGA_LIBRARY} == "MALBAC" ]]
then
        Adapter1="GTGAGTGATGGTTGAGGTAGTGTGGAG"
        Adapter2="CTCCACACTACCTCAACCATCACTCAC"
        adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2}"
elif [[ ${WGA_LIBRARY} == "PICOPLEX" ]]
then
        Adapter1="TGTGTTGGGTGTGTTTGG"
        Adapter2="CCAAACACACCCAACACA"
        Adapter3="TGTTGTGGGTTGTGTTGG"
        Adapter4="CCAACACAACCCACAACA"
        Adapter5="TGTGTTGGGTGTGTTTGG"
        Adapter6="CCAAACACACCCAACACA"
        adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2} -g Adapter5.2=${Adapter3} -G Adapter5.2=${Adapter3} -a Adapter3.2=${Adapter4} -A Adapter3.2=${Adapter4} -g Adapter5.3=${Adapter5} -G Adapter5.3=${Adapter5} -a Adapter3.3=${Adapter6} -A Adapter3.3=${Adapter6}" 
else 
        echo "CutAdapt not needed for ${WGA_LIBRARY}?"
        ln -s Processing/${sample}.trimmed_1.fastq.gz Processing/${sample}.trimmed2_1.fastq.gz
        ln -s Processing/${sample}.trimmed_2.fastq.gz Processing/${sample}.trimmed2_2.fastq.gz
        exit
fi

cutadapt \
	--minimum-length 70 \
    ${adapters_to_remove} \
    -o Processing/${sample}.trimmed2_1.fastq.gz \
    -p Processing/${sample}.trimmed2_2.fastq.gz \
    Processing/${sample}.trimmed_1.fastq.gz Processing/${sample}.trimmed_2.fastq.gz \
    > ./slurm_files/${sample}_CutadaptWGA.txt