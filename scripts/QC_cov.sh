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

for i in Processing/*.dedup.bam;do
j=$(basename "$i" .dedup.bam)

bedtools genomcov -ibam $i | grep "^genome" > Processing/$j.genome.bed

done
