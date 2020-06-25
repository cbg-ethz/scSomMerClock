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
        -n | --name)		shift
							name=$1
							;;
		-c | --chr)			shift
							chromosome=$1
							;;
        *)                  sample_bams+="$1 " 
    esac
    shift
done
echo ${sample_bams}

bams_in=$(echo ${sample_bams} | sed 's/ / -I /g')


echo ${sample_bams} | sed 's/ /\n/g' | sed 's/.recal.bam//g' \
| awk -v chr=${chromosome} '{print $0".recal.bam\t"$0".real."chr".bam"}' \
| sed 's/Processing\///' > ${name}.${chromosome}.map


java -Djava.io.tmpdir=Processing/ -Xmx25G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I ${bams_in} \
    -o ${name}.${chromosome}.intervals \
    -R ${REF} \
    -known ${INDELS} \
    -known ${INDELS2} \
    -L ${chromosome}

java -Djava.io.tmpdir=Processing/ -Xmx25G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known ${INDELS} \
    -known ${INDELS2} \
    -I ${bams_in} \
    -R ${REF} \
    -targetIntervals ${name}.${chromosome}.intervals \
    -L ${chromosome} \
    --nWayOut ${name}.${chromosome}.map \
    --maxReadsForRealignment 1000000