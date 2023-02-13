#!/bin/bash

tumor_in=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -gAD | --gnomAD )   shift
                            gnomAD=$1
                            ;;
        -n | --normal )     shift
                            normal_in="$1"
                            ;;
        *)                  tumor_in+="$1 "
    esac
    shift
done

set -Eeuxo pipefail

all_f1r2_input=`for chromosome in {1..22}; do
        printf -- "-I Calls/${chromosome}.f1r2.mutect.tar.gz "; done`
gatk LearnReadOrientationModel \
    ${all_f1r2_input} \
    --output Calls/read-orientation-model.tar.gz

## Normal bulk sample
normal_name=$(echo $normal_in | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" GetPileupSummaries \
    --input ${normal_in}  \
    --variant ${gnomAD} \
    --intervals ${gnomAD} \
    --output Calls/getpileupsummaries.${normal_name}.table

## Tumor bulk samples
for tumor in $tumor_in; do
    sample=$(echo $tumor | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
    gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" GetPileupSummaries \
        --input ${tumor}  \
        --variant ${gnomAD} \
        --intervals ${gnomAD} \
        --output Calls/getpileupsummaries.${sample}.table
    gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" CalculateContamination \
        --input Calls/getpileupsummaries.${sample}.table \
        --matched-normal Calls/getpileupsummaries.${normal_name}.table \
        --output Calls/${sample}.contamination.table
done