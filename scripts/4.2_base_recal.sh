#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Cell name
##$3: Reference genome file
cellname=$2
REF=$3

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" ApplyBQSR \
    -R ${REF} \
    -I Processing/${cellname}.dedup.bam \
    --bqsr Processing/${cellname}.recal.table \
    -O Processing/${cellname}.recal.bam