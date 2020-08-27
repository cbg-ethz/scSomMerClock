#!/bin/sh

module purge

minBaseQual=13
minDepth=10
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
        --mbq | --minBaseQual ) shift
                            minBaseQual=$1
                            ;;
        -md | --minDepth )  shift
                            minDepth=$1
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

cores=$(nproc)

python $SCcaller \
    --bam Processing/${cellname}.real.${chr}.bam \
    --fasta ${REF} \
    --output Calls/${cellname}.${chr}.sccaller.vcf \
    --snp_type dbsnp \
    --snp_in ${DBSNP} \
    --cpu_num ${cores} \
    --engine samtools \
    --bulk Processing/${bulk_normal}.real.${chr}.bam \
    --min_depth ${minDepth} \
    --mapq ${minBaseQual} \
    --minvar 2 
# && rm sc_${cellname}.real.${chr}.sccallerlab_01to-1.log
