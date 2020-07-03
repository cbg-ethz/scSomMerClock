#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Dataset name
##$3: Chromosome
##$4: Reference genome file
name=$2
chr=$3
REF=$4



samtools mpileup \
    -r ${chr} \
    -BQ0 \
    -d10000 \
    -f ${REF} \
    -q 30 \
    -b Processing/${chr}.bamspath.txt | \
monovar.py \
    -p 0.002 \
    -a 0.2 \
    -t 0.05 \
    -m 3 \
    -c 1 \
    -f ${REF} \
    -b Processing/${chr}.bamspath.txt \
    -o Calls/${name}.${chr}.monovar.vcf
