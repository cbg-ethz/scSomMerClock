#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user linkmonica@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 2
#SBATCH -t 60:00:00
#SBATCH --mem 50G
#SBATCH -p shared
#SBATCH --qos=shared

module load gcc/6.4.0
module load samtools
module load python/2.7.15
module load pysam
module load numpy/1.15.2-python-2.7.15


REF=/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa
DBSNP=/mnt/netapp1/posadalab/phylocancer/RESOURCES/dbsnp_138.b37.vcf
INDELS=/mnt/netapp1/posadalab/phylocancer/RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf

#SCcaller=/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Programs/SCcaller_2.0.0/sccaller_v2.0.0.py
SCcaller=/mnt/netapp1/posadalab/APPS/sccaller-2.0.0/SCcaller-2.0.0/sccaller_v2.0.0_lab.py

head -n $SLURM_ARRAY_TASK_ID $1 | tail -1 | while read -r cellnames;do

for chr in `seq 1 22` X Y;
do

python $SCcaller \
        --bam Processing/$cellnames.real.$chr.bam \
        --fasta $REF \
        --output Processing/$cellnames.real.$chr.sccallerlab.vcf \
        --snp_type dbsnp \
        --snp_in $DBSNP \
        --cpu_num 2 \
        --engine samtools \
        --bulk Processing/bulk.real.$chr.bam \
        --min_depth 1 \
        --minvar 0 \
        --mapq 30 \
	--bias 0.6 \
	--lamb 2000

done
done
