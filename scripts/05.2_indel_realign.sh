#!/bin/sh

module purge

bams_in=""
INDELS2=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -c | --chr)         shift
                            chr=$1
                            ;;
        -t | --target)      shift
                            target=$1
                            ;;
        -ma | --maps)        shift
                            maps=$1
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
    esac
    shift
done

set -Eeuxo pipefail
MEM_GB=$(awk '/Active:/ { printf "%.0f \n", $2/1024/1024 - 4}' /proc/meminfo)

java -Djava.io.tmpdir=Processing/ -Xmx${MEM_GB}G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known ${INDELS1} \
    ${INDELS2} \
    ${bams_in} \
    -R ${REF} \
    -targetIntervals ${target} \
    -L ${chr} \
    --nWayOut ${maps} \
    --maxReadsForRealignment 1000000
