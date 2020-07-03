#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Cell name
##$3: Reference genome file
##$4: DBSNP file
##$5: INDEL db file
cellname=$2
REF=$3
DBSNP=$4
INDELS=$5

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" BaseRecalibrator \
    -I Processing/${cellname}.dedup.bam \
    -O Processing/${cellname}.recal.table \
    -R ${REF} \
    --known-sites ${DBSNP} \
    --known-sites ${INDELS}