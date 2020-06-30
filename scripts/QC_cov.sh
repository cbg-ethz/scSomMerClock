#!/bin/sh

module purge
module load bedtools/2.28.0

cell=$1
bed=$2
SEQ=$3
WES_REF=$4

j=$(basename "$cell" .dedup.bam)
if [[ ${SEQ} == "WGA" ]]
then
    bedtools genomcov -ibam ${cell} | grep "^genome" > ${bed}
else
    bedtools coverage -a ${WES_REF} -b ${cell} -hist -sorted | grep "^all" > ${bed}
fi
