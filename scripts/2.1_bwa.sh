#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Sample name
##$3: Reference genome file
##$4: WGA Library string 
sample=$2
REF=$3
WGA_LIBRARY=$4

SM=$(echo ${sample} | cut -d "_" -f1)
PL=$(echo "ILLUMINA")
LB=$(echo "KAPA")
PU=`zcat Raw_Data/${sample}_1.fastq.gz | head -1 | sed 's/[:].*//' | sed 's/@//' \
    | sed 's/ /_/g'`
RG="@RG\\tID:${sample}\\tSM:${SM}\\tPL:${PL}\\tLB:${LB}\\tPU:${PU}"

if [[ ${WGA_LIBRARY} == "AMPLI-1" ]] || [[ ${WGA_LIBRARY} == "MALBAC" ]] || [[ ${WGA_LIBRARY} == "PICOPLEX" ]]
then    
    bwa mem -t 10 -R ${RG} ${REF} \
        Processing/${sample}.trimmed2_1.fastq.gz Processing/${sample}.trimmed2_2.fastq.gz \
        > Processing/${sample}.sam
else
    bwa mem -t 10 -R ${RG} ${REF} \
        Processing/${sample}.trimmed_1.fastq.gz Processing/${sample}.trimmed_2.fastq.gz \
        > Processing/${sample}.sam
fi