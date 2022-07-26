#!/bin/bash

sample_bams=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -o | --out)         shift
                            out_file=$1
                            ;;
        *)                  sample_bams+="$1 " 
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)

for sample in ${sample_bams}
do
    bcftools sort \
        --output-file ${sample}.gz \
        --output-type z \
        ${sample}
    bcftools index \
        --force \
        --threads ${cores} \
        ${sample}.gz
done

sorted_bams=$(echo "${sample_bams}" \
    | tr ' ' '\n' \
    | sort -V --field-separator=. --key=2 \
    | sed 's/vcf$/vcf\.gz/' \
    | tr '\n' ' '
)

bcftools concat \
    --output-type u \
    --threads ${cores} \
    ${sorted_bams} \
| bcftools filter \
    --exclude "TYPE!='snp' | FORMAT/SO!='True' | FILTER='multiple_genotype'" \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    -
bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}