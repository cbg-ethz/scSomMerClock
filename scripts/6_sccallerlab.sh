#!/bin/sh

module purge

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
    esac
    shift
done

[[ -z "$cellname" ]] && { echo "Error: Cellname not set"; exit 1; }
[[ -z "$chr" ]] && { echo "Error: Chromosome not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$SCcaller" ]] && { echo "Error: SCcaller exe not set"; exit 1; }
[[ -z "$bulk_normal" ]] && { echo "Error: Bulk normal not set"; exit 1; }
[[ -z "$DBSNP" ]] && { echo "Error: DBSNP file not set"; exit 1; }

python $SCcaller \
    --bam Processing/${cellname}.real.${chr}.bam \
    --fasta ${REF} \
    --output Calls/${cellname}.real.${chr}.sccallerlab.vcf \
    --snp_type dbsnp \
    --snp_in ${DBSNP} \
    --cpu_num 2 \
    --engine samtools \
    --bulk Processing/${bulk_normal}.real.${chr}.bam \
    --min_depth 1 \
    --minvar 0 \
    --mapq 30 \
    --bias 0.6 \
    --lamb 2000
