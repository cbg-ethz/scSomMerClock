#!/bin/sh

module purge

chr_vcfs=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)          shift
                                module load $1
                                ;;
        -o | --out)             shift
                                out_file=$1
                                ;;
        *)                      chr_vcfs+="$1 " 
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)

sorted_vcfs=$(echo ${chr_vcfs} \
    | tr ' ' '\n' \
    | sort -V \
    | tr '\n' ' '
)

bcftools concat \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    ${sorted_vcfs}
bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}
