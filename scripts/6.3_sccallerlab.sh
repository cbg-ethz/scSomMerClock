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
    bcftools query -l ${sample} \
        | awk '{print $0"\t"$0".sccaller" > "vcf_header.sccaller.tmp"}' \
    && bcftools reheader -s vcf_header.sccaller.tmp --threads ${cores} -o ${sample} ${sample}
done
rm vcf_header.sccaller.tmp

cores=$(nproc)
bcftools merge \
    --output ${out_file} \
    --output-type z \
    --merge both \
    --threads ${cores} \
    ${sample_bams}
