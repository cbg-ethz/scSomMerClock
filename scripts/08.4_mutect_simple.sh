#!/bin/sh

module purge

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
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)

stats=`for chromosome in {1..22}; do
    printf -- "--stats Calls/${chromosome}.mutect.vcf.stats "; done`
gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" MergeMutectStats \
        $stats \
        --output Calls/mutect.merged.stats \
&& gatk --java-options "-Xmx35G -Djava.io.tmpdir=Calls/" FilterMutectCalls \
        --reference ${REF} \
        --variant ${vcf_in} \
        --stats Calls/mutect.merged.stats \
        --create-output-variant-index \
        --output ${out_file}