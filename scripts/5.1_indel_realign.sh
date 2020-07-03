#!/bin/sh

##$1: Module names
module purge
module load $1
shift

sample_bams=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -n | --name)        shift
                            name=$1
                            ;;
        -c | --chr)         shift
                            chromosome=$1
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

bams_in=$(echo ${sample_bams} | sed 's/ / -I /g')

java -Djava.io.tmpdir=Processing/ -Xmx35G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I ${bams_in} \
    -o Reallignment/${name}.${chromosome}.intervals \
    -R ${REF} \
    -known ${INDELS1} \
    -known ${INDELS2} \
    -L ${chromosome}