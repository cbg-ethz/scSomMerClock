#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Cell name
##$3: Chromosome
##$4: Bulk normal "cell" name
##$5: Reference genome file
##$6: DBSNP file
##$7: SCCaller exe
cellnames=$2
chr=$3
bulk_normal=$4
REF=$5
DBSNP=$6
SCcaller=$7

python $SCcaller \
    --bam Processing/${cellnames}.real.${chr}.bam \
    --fasta ${REF} \
    --output Calls/${cellnames}.real.${chr}.sccallerlab.vcf \
    --snp_type dbsnp \
    --snp_in ${DBSNP} \
    --cpu_num 2 \
    --engine samtools \
    --bulk Processing/${bulk_normal}.real.${chr}.bam \
    --min_depth 1 \
    --minvar 0 \
    --mapq 30 \
    --bias 0.6 \
    --lamb 2000
