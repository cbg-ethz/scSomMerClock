#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user linkmonica@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 2
#SBATCH -t 100:00:00
#SBATCH --mem 50G
#SBATCH -p shared
#SBATCH --qos=shared

module purge
module load gatk/3.7-0-gcfedb67

REF=/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa
DBSNP=/mnt/netapp1/posadalab/phylocancer/RESOURCES/dbsnp_138.b37.vcf
INDELS=/mnt/netapp1/posadalab/phylocancer/RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf
INDELS2=/mnt/netapp1/posadalab/phylocancer/RESOURCES/1000G_phase1.indels.b37.vcf

sample_bams=$(ls /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/CRC09*.recal.bam)
bams_in=$(echo $sample_bams | sed 's/ / -I /g')

bulk_bams=$(ls /mnt/netapp1/posadalab/laurat/CRC09/workdir/1806_CRC09.XP_WGS_MG.1806KHX-0086.recal.bam)
bulk_in=$(echo $bulk_bams)

echo $sample_bams $bulk_bams |sed 's/ /\n/g' | sed 's/.recal.bam//g' | awk -v chr=$SLURM_ARRAY_TASK_ID '{print $0".recal.bam\t"$0".real."chr".bam"}' \
| sed 's/\/mnt\/lustre\/scratch\/home\/uvi\/be\/posadalustre\/CRC\/CRC09\///' \
| sed 's/\/mnt\/netapp1\/posadalab\/laurat\/CRC09\/workdir\///' \
| sed 's/\/mnt\/netapp1\/posadalab\/laurat\/CRC09\/workdir\//\/mnt\/lustre\/scratch\/home\/uvi\/be\/posadalustre\/CRC\/CRC09\//' > CRC09.$SLURM_ARRAY_TASK_ID.map


java -Djava.io.tmpdir=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/ -Xmx40G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-I $bams_in \
	-I $bulk_in \
	-o CRC09.$SLURM_ARRAY_TASK_ID.intervals \
	-R $REF \
	-known $INDELS \
	-known $INDELS2 \
	-L $SLURM_ARRAY_TASK_ID

java -Djava.io.tmpdir=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/ -Xmx40G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-known $INDELS \
	-known $INDELS2	\
	-I $bams_in \
	-I $bulk_in \
	-R $REF \
	-targetIntervals CRC09.$SLURM_ARRAY_TASK_ID.intervals \
	-L $SLURM_ARRAY_TASK_ID \
	--nWayOut CRC09.$SLURM_ARRAY_TASK_ID.map \
	--maxReadsForRealignment 1000000
