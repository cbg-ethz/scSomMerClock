#!/bin/sh

module purge
module load gcccore/6.4.0 cutadapt/1.18-python-3.7.0
module load gcc/6.4.0 bwa/0.7.17
module load picard/2.18.14

REF=/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa
##$1: file with names and files


head -n $SLURM_ARRAY_TASK_ID $1 | tail -1 | while read -r out f1;do

ID=${out}
SM=$(echo ${out} | cut -d "_" -f1)
PL=$(echo "ILLUMINA")
LB=$(echo "KAPA")
PU=`zcat /mnt/netapp2/posadalab2/DATA/${f1}_1.fastq.gz | head -1 | sed 's/[:].*//' | sed 's/@//'`
echo "SAMPLE: "${out}" ID: "${ID}" SM: "${SM}
RG="@RG\\tID:${ID}\\tSM:${SM}\\tPL:${PL}\\tLB:${LB}\\tPU:${PU}"

bwa mem -t 10 \
	-R ${RG} \
	$REF \
	/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/${out}.trimmed2_1.fastq.gz /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/${out}.trimmed2_2.fastq.gz \
	> /mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/${out}.sam

java -Xmx18g -jar $EBROOTPICARD/picard.jar SortSam \
	I=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/${out}.sam \
	TMP_DIR=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/ \
	O=/mnt/lustre/scratch/home/uvi/be/posadalustre/CRC/CRC09/${out}.sorted.bam \
	CREATE_INDEX=true \
	SORT_ORDER=coordinate

done
