#!/bin/sh

module purge
module load gatk/4.0.10.0

##$1: Cell name
##$2: Reference genome file
##$3: DBSNP file
##$4: INDEL db file
cellname=$1
REF=$2
DBSNP=$3
INDELS=$4

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" BaseRecalibrator \
    -I Processing/${cellname}.dedup.bam \
    -O Processing/${cellname}.recal.table \
    -R ${REF} \
    --known-sites ${DBSNP} \
    --known-sites ${INDELS}