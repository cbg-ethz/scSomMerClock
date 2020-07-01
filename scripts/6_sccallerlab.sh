#!/bin/sh

module purge
module load gcc/6.4.0
module load samtools
module load python/2.7.15
module load pysam
module load numpy/1.15.2-python-2.7.15

##$1: Cell name
##$2: Chromosome
##$3: Bulk normal "cell" name
##$4: Reference genome file
##$5: DBSNP file
##$6: SCCaller exe
cellnames=$1
chr=$2
bulk_normal=$3
REF=$4
DBSNP=$5
SCcaller=$6

python $SCcaller \
    --bam Processing/${cellnames}.real.${chr}.bam \
    --fasta ${REF} \
    --output Processing/${cellnames}.real.${chr}.sccallerlab.vcf \
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
