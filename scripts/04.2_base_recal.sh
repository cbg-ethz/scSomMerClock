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
        -i | --input )      shift
                            input=$1
                            ;;
        -o | --output )     shift
                            output=$1
                            ;;
        -t | --table )      shift
                            table=$1
                            ;;
    esac
    shift
done

set -Eeuxo pipefail

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" ApplyBQSR \
    -R ${REF} \
    -I ${input} \
    --bqsr ${table} \
    -O ${output}