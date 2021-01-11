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
MEM_GB=$(awk '/Active:/ { printf "%.0f \n", $2/1024/1024 - 2}' /proc/meminfo)

gatk --java-options "-Xmx${MEM_GB}G -Djava.io.tmpdir=Processing/" ApplyBQSR \
    -R ${REF} \
    -I ${input} \
    --bqsr ${table} \
    -O ${output}