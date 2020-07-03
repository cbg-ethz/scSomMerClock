#!/bin/sh

##$1: Module names
module purge
module load $1

##$2: Sample name
sample=$2

java -Xmx18g -jar $EBROOTPICARD/picard.jar SortSam \
    I=Processing/${sample}.sam \
    TMP_DIR=Processing/ \
    O=Processing/${sample}.sorted.bam \
    CREATE_INDEX=true \
    SORT_ORDER=coordinate