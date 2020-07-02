#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user linkmonica@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 3
#SBATCH -t 80:00:00
#SBATCH --mem 10G
#SBATCH -p shared
#SBATCH --qos=shared


module load gcc/6.4.0 monovar/140316-python-2.7.15 scipy/1.1.0-python-2.7.15

REF=$2

### use sbatch Monovar.sh filenames(bamspath.txt) Samplename

time(
ls Processing/*real.$SLURM_ARRAY_TASK_ID.bam | grep -v "Cancer" | grep -v "Normal" | grep -v "Polyps" > $1.$SLURM_ARRAY_TASK_ID.bamspath.txt 
samtools mpileup -r $SLURM_ARRAY_TASK_ID -BQ0 -d10000 -f $REF -q 30 -b $1.$SLURM_ARRAY_TASK_ID.bamspath.txt | \
monovar.py -p 0.002 -a 0.2 -t 0.05 -m 3 -c 1 -f $REF -b $1.$SLURM_ARRAY_TASK_ID.bamspath.txt \
-o Processing/$1.$SLURM_ARRAY_TASK_ID.monovar.vcf
)
