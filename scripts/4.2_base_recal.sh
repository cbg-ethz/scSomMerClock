#!/bin/sh

module purge
module load gatk/4.0.10.0

##$1: Cell name
##$2: Reference genome file
cellname=$1
REF=$2

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" ApplyBQSR \
    -R ${REF} \
    -I Processing/${cellname}.dedup.bam \
    --bqsr Processing/${cellname}.recal.table \
    -O Processing/${cellname}.recal.bam