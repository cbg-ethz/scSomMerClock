#!/bin/sh

module purge

sample_bams=""
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
        *)                  sample_bams+="$1 " 
    esac
    shift
done

[[ -z "$chr" ]] && { echo "Error: Chromosome not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$INDELS1" ]] && { echo "Error: First Indel file not set"; exit 1; }
[[ -z "$INDELS2" ]] && { echo "Error: Second Indel file not set"; exit 1; }

bams_in=$(echo ${sample_bams} | sed 's/ / -I /g')

java -Djava.io.tmpdir=Processing/ -Xmx35G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known ${INDELS1} \
    -known ${INDELS2} \
    -I ${bams_in} \
    -R ${REF} \
    -targetIntervals Reallignment/${chr}.intervals \
    -L ${chr} \
    --nWayOut Reallignment/${chr}.map \
    --maxReadsForRealignment 1000000
