#!/bin/sh

module purge

sample_bams=""
mutect_calls=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -o | --out)         shift
                            out_file=$1
                            ;;
        *)                  if [[ $1 == *"mutect"* ]]; then
                                mutect_calls=$1
                                sample_bams+="$1.tmp "
                            else
                                sample_bams+="$1 "
                            fi
    esac
    shift
done

[[ -z "$out_file" ]] && { echo "Error: Output file not set"; exit 1; }

cores=$(nproc)

if [ "$mutect_calls" != "" ]; then
    # Rename header column and index
    bcftools query -l $mutect_calls \
        | sed 's/\.mutect$//g' \
        | awk -F "[.]" '{print $0"\t"$1".mutect" > "vcf_header.mutect.tmp"}' \
    && bcftools reheader \
        --samples vcf_header.mutect.tmp \
        --threads ${cores} \
        ${mutect_calls} \
    | bcftools annotate \
        --remove FORMAT/AD \
        --output-type z \
        --output ${mutect_calls}.tmp \
        - \
    && bcftools index \
        --force \
        --threads ${cores} \
        ${mutect_calls}.tmp \
    && rm vcf_header.mutect.tmp
fi

bcftools merge \
    --output ${out_file} \
    --output-type z \
    --merge both \
    --threads ${cores} \
    ${sample_bams} \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}

if [ "$mutect_calls" != "" ]; then
    rm ${mutect_calls}.tmp 
fi