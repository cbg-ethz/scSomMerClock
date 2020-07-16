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

# && bgzip --stdout --index Calls/${cellname}.real.${chr}.sccaller.vcf \
    # > Calls/${cellname}.real.${chr}.sccaller.vcf.gz \

cores=$(nproc)

sorted_bams=$(echo "${sample_bams}" | sort -V --field-separator=. --key=3 | tr '\n' ' ') # | sed 's/$/.gz/'
bcftools concat \
    -output ${out_file} \
    -output-type z \
    --threads ${cores} \
    ${sorted_bams} \
&& bcftools index \
    --force \
    --threads ${cores} \
    ${out_file}

