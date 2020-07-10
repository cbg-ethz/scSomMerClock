#!/bin/sh

module purge

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -m | --module)      shift
                            module load $1
                            ;;
        --seq)              shift
                            SEQ=$1
                            ;;
        -o | --output )     shift
                            out_file=$1
                            ;;
        -i | --input )      shift
                            in_file=$1
                            ;;
        -e | --exome )      shift
                            WES_REF=$1
                            ;;
        -g | --genome )     shift
                            genome=$1
                            ;;                    
    esac
    shift
done

[[ -z "$in_file" ]] && { echo "Error: Input .bam file not set"; exit 1; }
[[ -z "$out_file" ]] && { echo "Error: Output .bed not set"; exit 1; }
[[ -z "$SEQ" ]] && { echo "Error: Sequencing platform not set"; exit 1; }
[[ -z "$WES_REF" && ${SEQ} == "WES" ]] && { echo "Error: Exome target file not set"; exit 1; }
[[ -z "$genome" && ${SEQ} == "WES" ]] && { echo "Error: Exome target file not set"; exit 1; }

if [[ ${SEQ} == "WGA" ]]
then
    bedtools genomcov -ibam ${in_file} | grep "^genome" > ${out_file}
else
    bedtools coverage -a ${WES_REF} -b ${in_file} -g ${genome} -hist -sorted \
        | grep "^all" > ${out_file}
fi
