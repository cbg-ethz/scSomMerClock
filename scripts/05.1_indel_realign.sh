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
                            ;;
    esac
    shift
done

set -Eeuxo pipefail

java -Djava.io.tmpdir=Processing/ -Xmx35G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    ${bams_in} \
    -o Realignment/${chr}.intervals \
    -R ${REF} \
    -known ${INDELS1} \
    -known ${INDELS2} \
    -L ${chr}