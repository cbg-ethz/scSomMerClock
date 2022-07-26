#!/bin/bash

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -i | --input)       shift
                            input=$1
                            ;;
        -o | --output)      shift
                            output=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -d | --dbsnp )      shift
                            DBSNP=$1
                            ;;
        -id | --indels )    shift
                            INDELS=$1
                            ;;
    esac
    shift
done

set -Eeuxo pipefail
MEM_GB=$(awk '/Active:/ { printf "%.0f \n", $2/1024/1024 - 2}' /proc/meminfo)

gatk --java-options "-Xmx${MEM_GB}G -Djava.io.tmpdir=Processing/" BaseRecalibrator \
    -I ${input} \
    -O ${output} \
    -R ${REF} \
    --known-sites ${DBSNP} \
    --known-sites ${INDELS}