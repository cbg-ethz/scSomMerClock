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

sorted_bams=$(echo "${sample_bams}"  | tr ' ' '\n' | sort -V | tr '\n' ' ') # | sed 's/$/.gz/'
echo $sorted_bams
bcftools concat \
    --output ${out_file}.tmp \
    --output-type z \
    --threads ${cores} \
    --no-version \
    ${sorted_bams}

# Rename header column and index
# bcftools query -l ${out_file}.tmp \
#     | sed 's/\.mutect$//g' \
#     | awk -F "[.]" '{print $0"\t"$1".mutect" > "vcf_header.mutect.tmp"}' \
# && bcftools reheader \
#     --samples vcf_header.mutect.tmp \
#     --threads ${cores} \
#     ${out_file}.tmp \
# |
bcftools annotate \
    --remove FORMAT/AD \
    --output-type z \
    --output ${out_file} \
    ${out_file}.tmp \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file} 
# && rm vcf_header.mutect.tmp ${out_file}.tmp



# java -cp $EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
#     -V $vcfs_in \
#     -R $REF \
#     -out Processing/$1.original.vcf 
