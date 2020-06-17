#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user linkmonica@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 2
#SBATCH -t 24:00:00
#SBATCH --mem 40G
#SBATCH -p amd-shared
#SBATCH --qos=amd-shared

module load bedtools/2.28.0
module load python/3.7.7
module load numpy/1.18.1-python-3.7.7

WES_REF=/home/uvi/be/nbo/data/data/resources/Agilent_SureSelect_Clinical_Research_Exome/S06588914_Covered.modified.bed
SEQ="WES"

if [[ ${SEQ} == "WGA" ]]
then
    for i in Processing/*.dedup.bam;do
        j=$(basename "$i" .dedup.bam)
        bedtools genomcov -ibam $i | grep "^genome" > Processing/$j.genome.bed
    done
else
    for i in Processing/*.dedup.bam;do
        j=$(basename "$i" .dedup.bam)
        bedtools intersect -a $i -b ${WES_REF} > Processing/$j.intersect.bam
        bedtools genomcov -ibam Processing/$j.intersect.bam | grep "^genome" > Processing/$j.genome.bed
    done
fi
