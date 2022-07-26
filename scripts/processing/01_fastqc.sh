#!/bin/bash

pair_end=true
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -s | --sample )     shift
                            sample=$1
                            ;;
        -l | --lib )        shift
                            WGA_LIBRARY=$1
                            ;;
        -se | --single-end )
                            pair_end=false
                            ;;
    esac
    shift
done

set -Eeuxo pipefail

if [ "${pair_end}" = true ]
then
    cutadapt \
    	--minimum-length 70 \
        --cores=0 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o Processing/${sample}.trimmed_1.fastq.gz \
        -p Processing/${sample}.trimmed_2.fastq.gz \
        Raw_Data/${sample}_1.fastq.gz Raw_Data/${sample}_2.fastq.gz
else
    cutadapt \
        --minimum-length 70 \
        --cores=0 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
        -o Processing/${sample}.trimmed_1.fastq.gz \
        Raw_Data/${sample}_1.fastq.gz
fi


if [[ ${WGA_LIBRARY} == "AMPLI-1" ]]
then
        Adapter1="GCTGTCAGTTAA"
        Adapter2="TTAACTGACAGCAGGAATCCCACT"
        if [ "${pair_end}" = true ]
        then
            adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2}"
        else
            adapters_to_remove="-g Adapter5=${Adapter1} -a Adapter3=${Adapter2}"
        fi
elif [[ ${WGA_LIBRARY} == "MALBAC" ]]
then
        Adapter1="GTGAGTGATGGTTGAGGTAGTGTGGAG"
        Adapter2="CTCCACACTACCTCAACCATCACTCAC"
        if [ "${pair_end}" = true ]
        then
            adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2}"
        else
            adapters_to_remove="-g Adapter5=${Adapter1} -a Adapter3=${Adapter2}"
        fi
elif [[ ${WGA_LIBRARY} == "PICOPLEX" ]]
then
        Adapter1="TGTGTTGGGTGTGTTTGG"
        Adapter2="CCAAACACACCCAACACA"
        Adapter3="TGTTGTGGGTTGTGTTGG"
        Adapter4="CCAACACAACCCACAACA"
        Adapter5="TGTGTTGGGTGTGTTTGG"
        Adapter6="CCAAACACACCCAACACA"
        if [ "${pair_end}" = true ]
        then
            adapters_to_remove="-g Adapter5=${Adapter1} -G Adapter5=${Adapter1} -a Adapter3=${Adapter2} -A Adapter3=${Adapter2} -g Adapter5.2=${Adapter3} -G Adapter5.2=${Adapter3} -a Adapter3.2=${Adapter4} -A Adapter3.2=${Adapter4} -g Adapter5.3=${Adapter5} -G Adapter5.3=${Adapter5} -a Adapter3.3=${Adapter6} -A Adapter3.3=${Adapter6}" 
        else
            adapters_to_remove="-g Adapter5=${Adapter1} -a Adapter3=${Adapter2} -g Adapter5.2=${Adapter3} -a Adapter3.2=${Adapter4} -g Adapter5.3=${Adapter5} -a Adapter3.3=${Adapter6}" 
        fi
else 
        echo "CutAdapt not needed for ${WGA_LIBRARY}?"
        exit 
fi


if [ "${pair_end}" = true ]
then
    cutadapt \
    	--minimum-length 70 \
        --cores=0 \
        ${adapters_to_remove} \
        -o Processing/${sample}.trimmed2_1.fastq.gz \
        -p Processing/${sample}.trimmed2_2.fastq.gz \
        Processing/${sample}.trimmed_1.fastq.gz \
        Processing/${sample}.trimmed_2.fastq.gz
else
    cutadapt \
        --minimum-length 70 \
        --cores=0 \
        ${adapters_to_remove} \
        -o Processing/${sample}.trimmed2_1.fastq.gz \
        Processing/${sample}.trimmed_1.fastq.gz
fi