#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        -s | --sample )     shift
                            cellname=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -d | --dbsnp )      shift
                            DBSNP=$1
                            ;;
        -i | --indels )   shift
                            INDELS=$1
                            ;;
    esac
    shift
done

[[ -z "$cellname" ]] && { echo "Error: Cellnames not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$INDELS" ]] && { echo "Error: Indel file not set"; exit 1; }
[[ -z "$DBSNP" ]] && { echo "Error: DBSNP file not set"; exit 1; }

gatk --java-options "-Xmx24G -Djava.io.tmpdir=Processing/" BaseRecalibrator \
    -I Processing/${cellname}.dedup.bam \
    -O Processing/${cellname}.recal.table \
    -R ${REF} \
    --known-sites ${DBSNP} \
    --known-sites ${INDELS}