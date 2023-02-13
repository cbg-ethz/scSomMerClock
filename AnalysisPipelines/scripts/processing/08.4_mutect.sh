#!/bin/bash

cont_tables=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
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

set -Eeuxo pipefail

stats=`for chromosome in {1..22}; do
    printf -- "--stats Calls/${chromosome}.mutect.vcf.stats "; done`

gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" MergeMutectStats \
    $stats \
    --output Calls/mutect.merged.stats
gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" FilterMutectCalls \
    --reference $REF \
    --variant ${vcf_in} \
    ${cont_tables} \
    --stats Calls/mutect.merged.stats \
    --ob-priors $rom \
    --create-output-variant-index \
    --output Calls/all.mutect.filtered.vcf.gz