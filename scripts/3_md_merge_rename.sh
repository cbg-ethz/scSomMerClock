#!/bin/sh

module purge
module load picard/2.18.14

##$1: Sample name
##$2: Corresponding cell name
f1=$1
f2=$2
cellnames=$3

if [ "$cellnames" == "" ]
then
    sample_bams=$(ls Processing/${f1}.sorted.bam)
    cellnames=$f2
else
    sample_bams=$(ls Processing/${f1}.sorted.bam Processing/${f2}.sorted.bam)
fi
bams_in=$(echo $sample_bams | sed 's/ / I=/g')
echo $bams_in
echo ${cellnames}

java -Xmx35g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=$bams_in \
    TMP_DIR=Processing/ \
    O=Processing/${cellnames}.dedup.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    M=Processing/Duplicates_${cellnames}.txt



