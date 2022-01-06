#!/bin/sh

module load samtools/1.9

for i in Processing/*.dedup.bam;do
    j=$(basename "$i" .dedup.bam)
    dep=$(samtools depth $i | awk '{sum+=$3} END { print "Average = ",sum/NR}' | sed 's/Average = //')
    echo $j $dep >> depth.txt
done
