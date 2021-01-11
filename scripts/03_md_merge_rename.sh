#!/bin/sh

module purge

bams_in=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -s | --sample )     shift
                            cellname=$1
                            ;;
        *)                  bams_in+="I=$1 " 
                            ;;
    esac
    shift
done

set -Eeuxo pipefail
MEM_GB=$(awk '/Active:/ { printf "%.0f \n", $2/1024/1024 - 2}' /proc/meminfo)

java -Xmx${MEM_GB}g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    ${bams_in} \
    TMP_DIR=Processing/ \
    O=Processing/${cellname}.dedup.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    M=Processing/Duplicates_${cellname}.txt



