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

for sample in ${sample_bams}
do
    bcftools view \
        --output-file ${sample}.gz \
        --output-type z \
        ${sample} \
    && bcftools index \
        --force \
        --threads ${cores} \
        ${sample}.gz
done

sorted_bams=$(echo "${sample_bams}" \
    | tr ' ' '\n' \
    | sort -V --field-separator=. --key=3 \
    | sed 's/vcf$/vcf\.gz/' \
    | tr '\n' ' ' \
)

bcftools concat \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    --no-version \
    ${sorted_bams} \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}

