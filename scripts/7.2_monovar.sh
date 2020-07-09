#!/bin/sh

module purge

sample_bams=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -i | --chr)         shift
                            chr=$1
                            ;;
        -o | --out)         shift
                            out_file=$1
                            ;;
        *)                  sample_bams+="$1 " 
    esac
    shift
done

[[ -z "$out_file" ]] && { echo "Error: Output file not set"; exit 1; }

cores=$(nproc)

# Rename header column (only sample, not chromosome), zip and index
sort -V sample_bams | while read -r f
do
    bcftools query -l ${f} | awk -F "[.]" '{print $0"\t"$1}' \
        | bcftools reheader -s - ${f} \
        && bcftools index -t ${f}
done

bcftools concat -o ${out_file} -O v --threads ${cores} ${sample_bams}