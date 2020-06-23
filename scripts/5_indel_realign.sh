#!/bin/sh

module purge
module load gatk/3.7-0-gcfedb67

sample_bams=""
while [ "$1" != "" ]; do
	key=$1
    case key in
        -r | --ref )        shift
                            REF=$1
                            ;;
        -i1 | --indels1 )   shift
                            INDELS=$1
                            ;;
        -i2 | --indels2 )   shift
                            INDEL2=$1
                            ;;
        *)                  sample_bams+="$1 " 
    esac
    shift
done
echo $sample_bams

sample_bams=$(ls Processing/*.recal.bam)
bams_in=$(echo $sample_bams | sed 's/ / -I /g')

echo $sample_bams | sed 's/ /\n/g' | sed 's/.recal.bam//g' \
| awk -v chr=$SLURM_ARRAY_TASK_ID '{print $0".recal.bam\t"$0".real."chr".bam"}' \
| sed 's/Processing\///' > $1.$SLURM_ARRAY_TASK_ID.map


java -Djava.io.tmpdir=Processing/ -Xmx25G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I $bams_in \
    -o $1.$SLURM_ARRAY_TASK_ID.intervals \
    -R $REF \
    -known $INDELS \
    -known $INDELS2 \
    -L $SLURM_ARRAY_TASK_ID

java -Djava.io.tmpdir=Processing/ -Xmx25G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known $INDELS \
    -known $INDELS2 \
    -I $bams_in \
    -R $REF \
    -targetIntervals $1.$SLURM_ARRAY_TASK_ID.intervals \
    -L $SLURM_ARRAY_TASK_ID \
    --nWayOut $1.$SLURM_ARRAY_TASK_ID.map \
    --maxReadsForRealignment 1000000