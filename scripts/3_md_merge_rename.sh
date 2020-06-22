#!/bin/sh

module purge
module load picard/2.18.14

REF=/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa
##$1: file with names and files

head -n $SLURM_ARRAY_TASK_ID $1 | tail -1 | while read -r cellnames;do
sample_bams=$(ls /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/${cellnames}-*.sorted.bam)
bams_in=$(echo $sample_bams | sed 's/ / I=/g')

java -Xmx35g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
	I=$bams_in \
	TMP_DIR=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/ \
	O=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/$cellnames.dedup.bam \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT \
	M=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/Duplicates_$cellnames.txt

done
