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

[[ -z "$cellname" ]] && { echo "Error: Cellname not set"; exit 1; }

java -Xmx32g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    ${bams_in} \
    TMP_DIR=Processing/ \
    O=Processing/${cellname}.dedup.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    M=Processing/Duplicates_${cellname}.txt


