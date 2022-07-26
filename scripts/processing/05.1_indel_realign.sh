#!/bin/bash

bams_in=""
INDELS2=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -o | --output)      shift
                            output=$1
                            ;;
        -c | --chr)         shift
                            chr=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -i1 | --indels1 )   shift
                            INDELS1=$1
                            ;;
        -i2 | --indels2 )   shift
                            INDELS2="-known $1 "
                            ;;
        *)                  bams_in+="-I $1 "
                            ;;
    esac
    shift
done

set -Eeuxo pipefail
MEM_GB=$(awk '/Active:/ { printf "%.0f \n", $2/1024/1024 - 4}' /proc/meminfo)
cores=$(nproc)

java -Djava.io.tmpdir=Processing/ -Xmx${MEM_GB}G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    ${bams_in} \
    -o ${output} \
    -R ${REF} \
    -known ${INDELS1} \
    ${INDELS2} \
    -L ${chr} \
    -nt ${cores}