#!/bin/sh

module purge

cont_tables=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -o | --out)         shift
                            out_file=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -i | --in)          shift
                            vcf_in=$1
                            ;;
        *)                  cont_tables+="--contamination-table $1 "
    esac
    shift
done

[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$out_file" ]] && { echo "Error: Output file not set"; exit 1; }
[[ -z "$vcf_in" ]] && { echo "Error: Input VCF file not set"; exit 1; }

cores=$(nproc)

gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" FilterMutectCalls \
        --reference ${REF} \
        --variant ${vcf_in} \
        --create-output-variant-index \
        --output ${out_file}