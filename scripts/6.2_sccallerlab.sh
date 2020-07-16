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

for sample in ${sample_bams}
do
    bcftools view \
        --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,21,X,Y \
        --output ${out_file}.gz \
        --output-type z
done

cores=$(nproc)

sorted_bams=$(echo "${sample_bams}" \
    | sort -V --field-separator=. --key=3 \
    | tr '\n' ' ' \
    | sed 's/$/.gz/' 
)
echo $sorted_bams
bcftools concat \
    --allow-overlaps \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    ${sorted_bams} \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}

