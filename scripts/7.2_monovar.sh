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
for sample in ${sample_bams}
do
    bcftools query -l ${f} | awk -F "[.]" '{print $0"\t"$1 > "vcf_header.monovar.tmp"}' \
        && bcftools reheader -s vcf_header.monovar.tmp --threads ${cores} -o ${f} ${f} \
        && bcftools index -t ${f}
done
rm vcf_header.monovar.tmp

sorted_bams=$(echo "${sample_bams}" | sort -V)
bcftools concat -o ${out_file} -O z --threads ${cores} ${sorted_bams}