#!/bin/sh

module purge

sample_bams=""
filter_str_in="'"
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)          shift
                                module load $1
                                ;;
        -fq | --filter_qual)    shift
                                filter_str_in+="QUAL>$1 &&"
                                ;;
        -fd | --filter_depth)   shift
                                filter_str_in+="FORMAT/DP>$1 &&"
                                ;;
        -o | --out)             shift
                                out_file=$1
                                ;;
        *)                      sample_bams+="$1 " 
    esac
    shift
done

filter_str="${filter_str_in%???}' "

[[ -z "$out_file" ]] && { echo "Error: Output file not set"; exit 1; }

cores=$(nproc)

# Rename header column (only sample, not chromosome), zip and index
for sample in ${sample_bams}
do
    chr=$(echo $sample | cut -d '/' -f 2 | cut -d '.' -f 1)
    bcftools query -l ${sample} \
        | sed 's/\.monovar$//g' \
        | awk -F "[.]" '{print $0"\t"$1".monovar" > "vcf_header.monovar.tmp"}' \
    && bcftools reheader \
        -s vcf_header.monovar.tmp \
        --threads ${cores} \
        ${sample} \
    | bcftools filter \
        --include ${filter_str} \
        --output ${sample}.tmp \
        --output-type z \
        --regions ${chr} \
        --threads ${cores} \
        - \
    && mv ${sample}.tmp ${sample} \
    && grep '^#\<contig' ${sample} \
        || sed -i "/^#CHROM.*/i ##contig=<ID=$chr,eta=-1>" ${sample}
done
rm vcf_header.monovar.tmp

sorted_bams=$(echo ${sample_bams} | tr ' ' '\n' | sort -V | tr '\n' ' ') # | sed 's/$/.gz/'
bcftools concat \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    ${sorted_bams} \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}
