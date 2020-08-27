#!/bin/sh

module purge

monovar=monovar.py
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
        -e | --exe )        shift
                            monovar=$1
                            ;; 
    esac
    shift
done

[[ -z "$chr" ]] && { echo "Error: Chromosome not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }

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
    -b Processing/${chr}.bamspath.txt \
    -f ${REF} \
    -o Calls/${chr}.monovar.vcf \
    -t 0.05 \
    -p 0.002 \
    -a 0.2 \
    -m ${cores} \
    -c 1