#!/bin/sh

module purge
module load gatk/4.0.10.0

REF=/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa
DBSNP=/mnt/netapp1/posadalab/phylocancer/RESOURCES/dbsnp_138.b37.vcf
INDELS=/mnt/netapp1/posadalab/phylocancer/RESOURCES/Mills_and_1000G_gold_standard.indels.b37.vcf

##$1: file with Cell/Sample names

head -n $SLURM_ARRAY_TASK_ID $1 | tail -1 | while read -r cellname;do

gatk --java-options "-Xmx18G -Djava.io.tmpdir=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/" BaseRecalibrator \
	-I /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/$cellname.dedup.bam \
	-O /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/$cellname.recal.table \
	-R $REF \
	--known-sites $DBSNP \
	--known-sites $INDELS


gatk --java-options "-Xmx18G -Djava.io.tmpdir=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/" ApplyBQSR \
	-R $REF \
	-I /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/$cellname.dedup.bam \
	--bqsr /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/$cellname.recal.table \
	-O /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/$cellname.recal.bam	

done
