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


for i in Processing/*.dedup.bam;do
    j=$(basename "$i" .dedup.bam)
    if [[ ${SEQ} == "WGA" ]]
    then
        bedtools genomcov -ibam $i | grep "^genome" > Processing/$j.genome.bed
    else
        bedtools coverage -a ${WES_REF} -b $i -hist > Processing/$j.genome.bed
    fi
done
