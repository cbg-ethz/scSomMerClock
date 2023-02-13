#!/bin/bash

bams_in=""
normal_in=""
while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -c | --chr)         shift
                            chr=$1
                            ;;
        -r | --ref )        shift
                            REF=$1
                            ;;
        -g | --germline )   shift
                            GERM_RES=$1
                            ;;
        -p | --pon )        shift
                            PON=$1
                            ;;
        -n | --normal )     shift
                            normal_in="-normal $1 "
                            ;;
        *)                  bams_in+="-I $1 "
    esac
    shift
done

set -Eeuxo pipefail

cores=$(nproc)

for bam in $bams_in; do
    if [[ "$bam" == "-I" ]]; then
        continue
    fi
    sample="$(cut -d '.' -f 1 <<< "${bam##*/}" )"
    java -Xmx28g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
        I=$bam \
        O=$bam.temp \
        RGID=$sample \
        RGSM=$sample \
        RGLB=KAPA \
        RGPL=ILLUMINA \
        RGPU=unit1
    mv $bam.temp $bam
    java -Xmx28g -jar $EBROOTPICARD/picard.jar BuildBamIndex \
        I=$bam
done

gatk Mutect2 \
    --reference ${REF} \
    ${bams_in} \
    ${normal_in} \
    --intervals ${chr} \
    --panel-of-normals ${PON} \
    --germline-resource ${GERM_RES} \
    --af-of-alleles-not-in-resource 0.0000025 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --output Calls/${chr}.mutect.vcf \
    --bam-output Calls/${chr}.mutect.bam\
    --f1r2-tar-gz Calls/${chr}.f1r2.mutect.tar.gz \
    --native-pair-hmm-threads ${cores}
    # --genotype-pon-sites true 