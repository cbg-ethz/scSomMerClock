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

bcftools merge \
    --output ${out_file}.tmp \
    --output-type z \
    --merge both \
    --threads ${cores} \
    ${sample_bams}

bcftools query -l ${out_file}.tmp \
    | sed 's/\.sccaller$//g' \
    | awk '{print $0"\t"$0".sccaller" > "vcf_header.sccaller.tmp"}' \
&& bcftools reheader \
    --samples vcf_header.sccaller.tmp \
    --threads ${cores} \
    --output ${out_file} \
    ${out_file}.tmp \
&& rm vcf_header.sccaller.tmp ${out_file}.tmp