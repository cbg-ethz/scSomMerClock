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
        -rom | --read-orientation-model ) shift
                            rom=$1
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
[[ -z "$rom" ]] && { echo "Error: Read-Orientation-Model not set"; exit 1; }

cores=$(nproc)

stats=`for chromosome in {1..22}; do
    printf -- "--stats Calls/${chromosome}.mutect.vcf.stats "; done`

gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" MergeMutectStats \
        $stats \
        --output Calls/mutect.merged.stats \
&& gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" FilterMutectCalls \
        --reference $REF \
        --variant ${vcf_in} \
        ${cont_tables} \
        --stats Calls/mutect.merged.stats \
        --ob-priors $rom \
        --create-output-variant-index \
        --output Calls/all.mutect.filtered.vcf