#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -s | --sample )     shift
                            cellname=$1
                            ;;
    esac
    shift
done

[[ -z "$cellname" ]] && { echo "Error: Cellname not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" ApplyBQSR \
    -R ${REF} \
    -I Processing/${cellname}.dedup.bam \
    --bqsr Processing/${cellname}.recal.table \
    -O Processing/${cellname}.recal.bam