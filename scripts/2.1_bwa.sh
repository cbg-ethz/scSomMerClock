#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -s | --sample )     shift
                            sample=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -l | --lib )        shift
                            WGA_LIBRARY=$1
                            ;;
    esac
    shift
done

[[ -z "$sample" ]] && { echo "Error: Sample not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$WGA_LIBRARY" ]] && { echo "Error: WGA library not set"; exit 1; }


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