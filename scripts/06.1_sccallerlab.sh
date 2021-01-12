#!/bin/sh
module purge

minBaseQual=13
minDepth=10
bulk_normal=""
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
        -s | --sample )     shift
                            cellname=$1
                            ;;
        -e | --exe )        shift
                            SCcaller=$1
                            ;;
        -b | --bulk )       shift
                            bulk_normal=$1
                            ;;
        -d | --dbsnp )      shift
                            DBSNP=$1
                            ;;
        -md | --minDepth )  shift
                            minDepth=$1
                            ;;
    esac
    shift
done

set -Eeuxo pipefail

if [ -z ${bulk_normal} ]
    bulk=""
else
    bulk="--bulk Processing/${bulk_normal}.recal.${chr}.bam"
fi

cores=$(nproc)

python $SCcaller \
    --bam Processing/${cellname}.recal.${chr}.bam \
    --fasta ${REF} \
    --snp_type dbsnp \
    --snp_in ${DBSNP} \
    ${bulk} \
    --output Calls/${cellname}.${chr}.sccaller.vcf \
    --cpu_num ${cores} \
    --engine pysam \
    --min_depth ${minDepth} \
    --minvar 2 \
    --bulk_min_depth 0
# && rm sc_${cellname}.recal.${chr}.sccallerlab_01to-1.log
