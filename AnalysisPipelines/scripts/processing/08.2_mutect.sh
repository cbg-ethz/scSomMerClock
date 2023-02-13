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
sorted_bams=$(echo "${sample_bams}"  | tr ' ' '\n' | sort -V | tr '\n' ' ') # | sed 's/$/.gz/'

bcftools concat \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    --no-version \
    ${sorted_bams}
bcftools index \
    --tbi \
    --force \
    --threads ${cores} \
    ${out_file}