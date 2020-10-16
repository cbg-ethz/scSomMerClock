#!/bin/sh

module purge

bams_in=""
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
                            INDELS2=$1
                            ;;
        *)                  bams_in+="-I $1 "
    esac
    shift
done

set -Eeuxo pipefail

java -Djava.io.tmpdir=Processing/ -Xmx35G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known ${INDELS1} \
    -known ${INDELS2} \
    ${bams_in} \
    -R ${REF} \
    -targetIntervals ${target} \
    -L ${chr} \
    --nWayOut ${maps} \
    --maxReadsForRealignment 1000000
