#!/bin/sh

sample_bams=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -o | --out)             shift
                                out_file=$1
                                ;;
        *)                      sample_bams+="$1 " 
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)
sample_bams_ordered=$(echo ${sample_bams} | sort -V | tr '\n' ' ')

bcftools merge \
    --output ${out_file}.tmp \
    --output-type z \
    --merge both \
    --apply-filters .,PASS \
    --info-rules NS:sum \
    --threads ${cores} \
    ${sample_bams_ordered}

bcftools query -l ${out_file}.tmp \
    | sort -V \
    | sed 's/\.sccaller$//g' \
    | awk '{print $0"\t"$0".sccaller" > "vcf_header.sccaller.tmp"}'

bcftools reheader \
    --samples vcf_header.sccaller.tmp \
    --threads ${cores} \
    --output ${out_file} \
    ${out_file}.tmp

bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}
    
rm vcf_header.sccaller.tmp ${out_file}.tmp