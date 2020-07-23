#!/bin/sh

module purge

sample_bams=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)          shift
                                module load $1
                                ;;
        -o | --out)             shift
                                out_file=$1
                                ;;
        *)                      sample_bams+="$1 " 
    esac
    shift
done

[[ -z "$out_file" ]] && { echo "Error: Output file not set"; exit 1; }

cores=$(nproc)

# Rename header column (only sample, not chromosome), zip and index
for sample in ${sample_bams}
do
    chr=$(echo $sample | cut -d '/' -f 2 | cut -d '.' -f 1)
    sample_order=$(bcftools query -l ${sample} \
        | sort -V \
        | sed 's/\.monovar$//g' \
        | tr '\n' ','
    )
    echo ${sample_order%?} \
        | tr ',' '\n' \
        | awk -F "[.]" '{print $0"\t"$1".monovar" > "vcf_header.monovar.tmp"}'
    
    bcftools view \
        --samples ${sample_order%?} \
        --output-type z \
        ${sample} \
    | bcftools reheader \
        --samples vcf_header.monovar.tmp \
        --threads ${cores} \
        --output ${sample}.gz \
        - \
    && bcftools index \
        --force \
        --threads ${cores} \
        ${sample}.gz
    && grep '^#\<contig' ${sample}.gz \
        || sed -i "/^#CHROM.*/i ##contig=<ID=$chr,eta=-1>" ${sample}.gz
done
rm vcf_header.monovar.tmp

sorted_bams=$(echo ${sample_bams} \
    | tr ' ' '\n' \
    | sort -V \
    | sed 's/$/.gz/' \
    | tr '\n' ' '
)

bcftools merge \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    ${sorted_bams} \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}
