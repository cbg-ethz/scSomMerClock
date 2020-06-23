#!/bin/sh

module purge
module load gcccore/6.4.0 cutadapt/1.18-python-3.7.0
module load gcc/6.4.0 bwa/0.7.17
module load picard/2.18.14

##$1: Sample name
##$2: Reference genome file
out=$1
REF=$2

ID=${out}
SM=$(echo ${out} | cut -d "_" -f1)
PL=$(echo "ILLUMINA")
LB=$(echo "KAPA")
PU=`zcat Raw_Data/${out}_1.fastq.gz | head -1 | sed 's/[:].*//' | sed 's/@//' | sed 's/ /_/g'`
echo "SAMPLE: "${out}" ID: "${ID}" SM: "${SM}
RG="@RG\\tID:${ID}\\tSM:${SM}\\tPL:${PL}\\tLB:${LB}\\tPU:${PU}"
echo ${RG}

bwa mem -t 10 \
    -R ${RG} \
    $REF \
    Processing/${out}.trimmed_1.fastq.gz Processing/${out}.trimmed_2.fastq.gz \
    > Processing/${out}.sam

java -Xmx18g -jar $EBROOTPICARD/picard.jar SortSam \
        I=Processing/${out}.sam \
        TMP_DIR=Processing/ \
        O=Processing/${out}.sorted.bam \
        CREATE_INDEX=true \
        SORT_ORDER=coordinate