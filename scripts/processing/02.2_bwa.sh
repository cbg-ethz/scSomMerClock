#!/bin/sh

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -s | --sample )     shift
                            sample=$1
                            ;;
    esac
    shift
done

set -Eeuxo pipefail
MEM_GB=$(awk '/Active:/ { printf "%.0f \n", $2/1024/1024 - 2}' /proc/meminfo)

java -Xmx${MEM_GB}g -jar $EBROOTPICARD/picard.jar SortSam \
    I=Processing/${sample}.sam \
    TMP_DIR=Processing/ \
    O=Processing/${sample}.sorted.bam \
    CREATE_INDEX=true \
    SORT_ORDER=coordinate

java -Xmx${MEM_GB}g -jar $EBROOTPICARD/picard.jar ValidateSamFile \
    I=Processing/${sample}.sorted.bam \
    IGNORE_WARNINGS=true \
    MODE=VERBOSE \
    O=sanity_check_alignment.txt