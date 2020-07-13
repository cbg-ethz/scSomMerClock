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
    esac
    shift
done

[[ -z "$chr" ]] && { echo "Error: Chromosome not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }

cores=$(nproc)

samtools mpileup \
    -r ${chr} \
    -BQ0 \
    -d10000 \
    -f ${REF} \
    -q 30 \
    -b Processing/${chr}.bamspath.txt \
| monovar.py \
    -p 0.002 \
    -a 0.2 \
    -t 0.05 \
    -m ${cores} \
    -c 1 \
    -f ${REF} \
    -b Processing/${chr}.bamspath.txt \
    -o Calls/${chr}.monovar.vcf \
&& bgzip -c -i Calls/${chr}.monovar.vcf > Calls/${chr}.monovar.vcf.gz