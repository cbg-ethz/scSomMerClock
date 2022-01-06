#!/bin/sh

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
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

set -Eeuxo pipefail

if [[ ${SEQ} == "WGA" ]]
then
    bedtools genomcov -ibam ${in_file} | grep "^genome" > ${out_file}
else
    bedtools coverage -a ${WES_REF} -b ${in_file} -g ${genome} -hist -sorted \
        | grep "^all" > ${out_file}
fi
