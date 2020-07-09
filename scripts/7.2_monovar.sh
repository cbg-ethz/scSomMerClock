#!/bin/sh

module purge

sample_bams=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
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
for f in ${sample_bams}
do
    bcftools query -l ${f} | awk -F "[.]" '{print $0"\t"$1}' \
        | bcftools reheader -s - ${f} \
        && bcftools index -t ${f}
done

sorted_bams=$(echo "${sample_bams}" | sort -V)
bcftools concat -o ${out_file} -O v --threads ${cores} ${sorted_bams}