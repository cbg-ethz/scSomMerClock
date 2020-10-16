#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -i | --input)       shift
                            input=$1
                            ;;
        -o | --output)      shift
                            output=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -d | --dbsnp )      shift
                            DBSNP=$1
                            ;;
        -id | --indels )    shift
                            INDELS=$1
                            ;;
    esac
    shift
done

set -Eeuxo pipefail

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" BaseRecalibrator \
    -I ${input} \
    -O ${output} \
    -R ${REF} \
    --known-sites ${DBSNP} \
    --known-sites ${INDELS}