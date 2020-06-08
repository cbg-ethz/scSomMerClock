#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user linkmonica@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 2
#SBATCH -t 10:00:00
#SBATCH --mem 10G
#SBATCH -p amd-shared
#SBATCH --qos=amd-shared

module purge
module load gcccore/6.4.0 cutadapt/1.18-python-3.7.0

REF=/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa

WGA_LIBRARY="AMPLI-1"

if [[ ${WGA_LIBRARY} == "AMPLI-1" ]]
then
	Adapter1="GCTGTCAGTTAA"
        Adapter2="TTAACTGACAGCAGGAATCCCACT"
	adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2}"
elif [[ ${WGA_LIBRARY} == "MALBAC" ]]
then
        Adapter1="GTGAGTGATGGTTGAGGTAGTGTGGAG"
	Adapter2="CTCCACACTACCTCAACCATCACTCAC"
	adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2}"
elif [[ ${WGA_LIBRARY} == "PICOPLEX" ]]
then
        Adapter1="TGTGTTGGGTGTGTTTGG"
        Adapter2="CCAAACACACCCAACACA"
	Adapter3="TGTTGTGGGTTGTGTTGG"
	Adapter4="CCAACACAACCCACAACA"
	Adapter5="TGTGTTGGGTGTGTTTGG"
	Adapter6="CCAAACACACCCAACACA"
	adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2} -g Adapter5.2=${Adapter3} -G Adapter5.2=${Adapter3} -a Adapter3.2=${Adapter4} -A Adapter3.2=${Adapter4} -g Adapter5.3=${Adapter5} -G Adapter5.3=${Adapter5} -a Adapter3.3=${Adapter6} -A Adapter3.3=${Adapter6}" 
else 
	echo "CutAdapt not needed for ${WGA_LIBRARY}?"
	ln -s ${WORKDIR}/${SAMPLE}.trimmed_1.fastq.gz ${WORKDIR}/${SAMPLE}.trimmed2_1.fastq.gz
	ln -s ${WORKDIR}/${SAMPLE}.trimmed_2.fastq.gz ${WORKDIR}/${SAMPLE}.trimmed2_2.fastq.gz
	exit
fi


head -n $SLURM_ARRAY_TASK_ID $1 | tail -1 | while read -r fastq;do
cutadapt --minimum-length 70 \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-o Processing/${fastq}.trimmed_1.fastq.gz -p Processing/${fastq}.trimmed_2.fastq.gz \
	Raw_Data/${fastq}_1.fastq.gz Raw_Data/${fastq}_2.fastq.gz > ./slurm_files/${fastq}_Cutadapt.txt

cutadapt --minimum-length 70 \
	${adapters_to_remove} \
	-o Processing/${fastq}.trimmed2_1.fastq.gz -p Processing/${fastq}.trimmed2_2.fastq.gz \
	Processing/${fastq}.trimmed_1.fastq.gz Processing/${fastq}.trimmed_2.fastq.gz > ./slurm_files/${fastq}_CutadaptWGA.txt

done
