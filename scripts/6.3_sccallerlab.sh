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
                                filter_str_in+="QUAL<$1 &&"
                                ;;
        -fd | --filter_depth)   shift
                                filter_str_in+="FORMAT/AD[:0-1]<$1 &&"
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

bcftools merge \
    --output ${out_file}.tmp \
    --output-type z \
    --merge both \
    --threads ${cores} \
    ${sample_bams}

bcftools query -l ${out_file}.tmp \
    | sort -V \
    | sed 's/\.sccaller$//g' \
    | awk '{print $0"\t"$0".sccaller" > "vcf_header.sccaller.tmp"}' \
&& bcftools reheader \
    --samples vcf_header.sccaller.tmp \
    --threads ${cores} \
    ${out_file}.tmp \
| bcftools filter \
    --exclude ${filter_str} \
    --output ${out_file} \
    --output-type z \
    --threads ${cores} \
    - \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file} 
# && rm vcf_header.sccaller.tmp ${out_file}.tmp