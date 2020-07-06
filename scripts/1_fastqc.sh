#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -s | --sample )     shift
                            sample=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -l | --lib )        shift
                            WGA_LIBRARY=$1
                            ;;
    esac
    shift
done

[[ -z "$sample" ]] && { echo "Error: Sample not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$WGA_LIBRARY" ]] && { echo "Error: WGA library not set"; exit 1; }


mkdir -p Adapter_Cutting
cutadapt \
	--minimum-length 70 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    -o Processing/${sample}.trimmed_1.fastq.gz \
    -p Processing/${sample}.trimmed_2.fastq.gz \
    Raw_Data/${sample}_1.fastq.gz Raw_Data/${sample}_2.fastq.gz \
    > Adapter_Cutting/${sample}_Cutadapt.txt


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
        ln -s Processing/${sample}.trimmed_1.fastq.gz Processing/${sample}.trimmed2_1.fastq.gz
        ln -s Processing/${sample}.trimmed_2.fastq.gz Processing/${sample}.trimmed2_2.fastq.gz
        exit
fi

cutadapt \
	--minimum-length 70 \
    ${adapters_to_remove} \
    -o Processing/${sample}.trimmed2_1.fastq.gz \
    -p Processing/${sample}.trimmed2_2.fastq.gz \
    Processing/${sample}.trimmed_1.fastq.gz Processing/${sample}.trimmed_2.fastq.gz \
    > Adapter_Cutting/${sample}_CutadaptWGA.txt