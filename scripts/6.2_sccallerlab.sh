#!/bin/sh

module purge

sample_bams=""
min_depth=10
min_qual=30
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -o | --out)         shift
                            out_file=$1
                            ;;
        -md | --min-depth)  shift
                            min_depth=$1
                            ;;
        -mq | --min-qual)   shift
                            min_qual=$1
                            ;;
        *)                  sample_bams+="$1 " 
    esac
    shift
done

[[ -z "$out_file" ]] && { echo "Error: Output file not set"; exit 1; }

filter_str='FORMAT/SO=="True"'" & QUAL >= ${min_qual} & FORMAT/AD[0:0] + FORMAT/AD[0:1] >= ${min_depth}"
cores=$(nproc)

for sample in ${sample_bams}
do
    bcftools filter \
        --include "${filter_str}" \
        --threads ${cores} \
        ${sample} \
    | bcftools sort  \
        --output-file ${sample}.gz \
        --output-type z \
        - \
    && bcftools index \
        --force \
        --threads ${cores} \
        ${sample}.gz
done

sorted_bams=$(echo "${sample_bams}" \
    | tr ' ' '\n' \
    | sort -V --field-separator=. --key=2 \
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

