#!/bin/sh

module purge

monovar=monovar.py
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
        -e | --exe )        shift
                            monovar=$1
                            ;; 
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)

samtools mpileup \
    --region ${chr} \
    --no-BAQ \
    --min-BQ 13 \
    --max-depth 10000 \
    --fasta-ref ${REF} \
    --min-MQ 40 \
    --bam-list Processing/${chr}.bamspath.txt \
    | ${monovar} \
        --bam_file_list Processing/${chr}.bamspath.txt \
        --ref_file ${REF} \
        --output Calls/${chr}.monovar.vcf \
        --threshold 0.05 \
        --pe 0.002 \
        --pad 0.2 \
        --cpus ${cores} \
        --CF_flag 0