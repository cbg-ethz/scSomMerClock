#!/bin/bash


sample_bams=""
mutect_calls=""
out_dir='.'
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -o | --out)         shift
                            out_dir=$1
                            ;;
        *)                  if [[ $1 == *"mutect"* ]]; then
                                mutect_calls=$1
                                sample_bams+="$1.tmp.gz "
                            else
                                sample_bams+="$1 "
                            fi
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)

if [ "$mutect_calls" != "" ]; then
    # Rename header column and index
    bcftools query -l $mutect_calls \
        | sed 's/\.mutect$//g' \
        | awk -F "[.]" '{print $0"\t"$1".mutect" > "vcf_header.mutect.tmp"}'
    bcftools reheader \
        --samples vcf_header.mutect.tmp \
        --threads ${cores} \
        ${mutect_calls} \
        | bcftools annotate \
            --output-type v \
            --output ${mutect_calls}.tmp \
            -
    rm vcf_header.mutect.tmp
    sed -i -e 's/ID=AD,Number=R,Type/ID=AD,Number=.,Type/g' ${mutect_calls}.tmp

    bcftools view ${mutect_calls}.tmp -O z -o ${mutect_calls}.tmp.gz
    bcftools index \
        --force \
        --threads ${cores} \
        ${mutect_calls}.tmp.gz
fi

bcftools merge \
    --output-type u \
    --merge both \
    --threads ${cores} \
    ${sample_bams} \
    | bcftools filter \
        --include 'N_PASS(GT="alt") > 0' \
        --threads ${cores} \
        --output-type z \
        --output ${out_dir}/all.vcf.gz \
        - 
bcftools index \
    --force \
    --threads ${cores} \
    ${out_dir}/all.vcf.gz

if [ "$mutect_calls" != "" ]; then
    rm ${mutect_calls}.tmp ${mutect_calls}.tmp.gz
fi

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 'X' 'Y'; do 
    bcftools view \
        --regions ${chr} \
        --types snps \
        --output-file ${out_dir}/all.${chr}.vcf.gz \
        --output-type z \
        ${out_dir}/all.vcf.gz
    bcftools index \
        --force \
        --threads ${cores} \
        ${out_dir}/all.${chr}.vcf.gz
done