#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Sorted bam files of Sample [1 or 2]
##$3: Corresponding cell name
f1=$2
f2=$3
cellnames=$4

if [ "${cellnames}" == "" ]
then
    sample_bams=${f1}
    cellnames=${f2}
else
    sample_bams="${f1} ${f2}"
fi
bams_in=$(echo ${sample_bams} | sed 's/ / I=/g')

java -Xmx35g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=${bams_in} \
    TMP_DIR=Processing/ \
    O=Processing/${cellnames}.dedup.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    M=Processing/Duplicates_${cellnames}.txt



