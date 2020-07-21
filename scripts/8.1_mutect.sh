#!/bin/sh

module purge

bams_in=""
normal_in=""
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
        -g | --germline )   shift
                            GERM_RES=$1
                            ;;
        -p | --pon )        shift
                            PON=$1
                            ;;
        -n | --normal )     shift
                            normal_in="-normal $1 "
                            ;;
         *)                 bams_in+="-I $1 "
    esac
    shift
done

[[ -z "$chr" ]] && { echo "Error: Chromosome not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$GERM_RES" ]] && { echo "Error: Germline resources not set"; exit 1; }
[[ -z "$PON" ]] && { echo "Error: Panel-Of-Normals not set"; exit 1; }

cores=$(nproc)

gatk Mutect2 \
    --reference ${REF} \
    ${bams_in} \
    ${normal_in} \
    --intervals ${chr} \
    --genotype-pon-sites true \
    --germline-resource ${GERM_RES} \
    --panel-of-normals ${PON} \
    --output Calls/${chr}.unfiltered.mutect.vcf \
    --f1r2-tar-gz Calls/${chr}.f1r2.mutect.tar.gz \
    --native-pair-hmm-threads ${cores} \
&& gatk FilterMutectCalls \
    --variant Calls/${chr}.unfiltered.mutect.vcf \
    --output Calls/${chr}.filtered.mutect.vcf \
    --reference ${REF}


# gatk LearnReadOrientationModel \
#     --input Calls/${chr}.f1r2.mutect.tar.gz \
#     --output Calls/${chr}.rom.mutect.tar.gz \
# && gatk GetPileupSummaries \
#     --input Calls/${chr}.rom.mutect.tar.gz  \
#     --variant Calls/${chr}.unfiltered.mutect.vcf \
#     --intervals ${chr} \
#     --output Calls/${chr}.getpileupsummaries.table \
# && gatk CalculateContamination \
#     --input Calls/${chr}.getpileupsummaries.table \
#     --matched-normal normal-pileups.table \
#     --output Calls/${chr}.calculatecontamination.table \
# && gatk FilterMutectCalls \
#     --variant Calls/${chr}.unfiltered.mutect.vcf \
#     --contamination-table Calls/${chr}.calculatecontamination.table \
#     --ob-priors Calls/${chr}.rom.mutect.tar.gz \
#     --output Calls/${chr}.filtered.mutect.vcf \
#     --reference ${REF}
