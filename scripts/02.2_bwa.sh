#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -s | --sample )     shift
                            sample=$1
                            ;;
    esac
    shift
done

[[ -z "$sample" ]] && { echo "Error: Sample not set"; exit 1; }

java -Xmx32g -jar $EBROOTPICARD/picard.jar SortSam \
    I=Processing/${sample}.sam \
    TMP_DIR=Processing/ \
    O=Processing/${sample}.sorted.bam \
    CREATE_INDEX=true \
    SORT_ORDER=coordinate

java -Xmx32g -jar $EBROOTPICARD/picard.jar ValidateSamFile \
    I=Processing/${sample}.sorted.bam \
    IGNORE_WARNINGS=true \
    MODE=VERBOSE \
    O=sanity_check_alignment.txt