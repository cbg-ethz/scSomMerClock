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
    --output Calls/${chr}.mutect.vcf \
    --f1r2-tar-gz Calls/${chr}.f1r2.mutect.tar.gz \
    --native-pair-hmm-threads ${cores} \
&& gatk FilterMutectCalls \
    --variant Calls/${chr}.mutect.vcf \
    --output Calls/${chr}.filtered.mutect.vcf \
    --reference ${REF}


# gatk LearnReadOrientationModel $all_f1r2_input -O Processing/read-orientation-model.tar.gz


# && gatk GetPileupSummaries \
#     --input Calls/${chr}.rom.mutect.tar.gz  \
#     --variant Calls/${chr}.mutect.vcf \
#     --intervals ${chr} \
#     --output Calls/${chr}.getpileupsummaries.table \
# && gatk CalculateContamination \
#     --input Calls/${chr}.getpileupsummaries.table \
#     --matched-normal normal-pileups.table \
#     --output Calls/${chr}.calculatecontamination.table \
# && gatk FilterMutectCalls \
#     --variant Calls/${chr}.mutect.vcf \
#     --contamination-table Calls/${chr}.calculatecontamination.table \
#     --ob-priors Calls/${chr}.rom.mutect.tar.gz \
#     --output Calls/${chr}.filtered.mutect.vcf \
#     --reference ${REF}

# ------------------------------------------------------------------------------

# module load gatk/4.1.1.0


# chr_vcfs=$(ls Processing/Wu61.*.mutect.vcf)
# vcfs_in=$(echo $chr_vcfs | sed 's/ / -V /g')


# module purge
# module load gatk/3.7-0-gcfedb67

# java -cp $EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
#         -V $vcfs_in \
#         -R $REF \
#         -out Processing/Wu61.original.vcf 

# module purge
# module load gatk/4.1.1.0

# stats_files=$(ls Processing/Wu61.*.mutect.vcf.stats)
# stats=$(echo $stats_files | sed 's/ / -stats /g')

# gatk --java-options "-Xmx28G -Djava.io.tmpdir=Processing/" MergeMutectStats \
#         -stats $stats \
#         -O Processing/Wu61.merged.stats

# gatk --java-options "-Xmx28G -Djava.io.tmpdir=Processing/" FilterMutectCalls \
#         -R $REF \
#         -V Processing//Wu61.original.vcf \
#         --contamination-table Processing/CRC0907-Cancer.contamination.table \
#         --contamination-table Processing/CRC0907-polyps.contamination.table \
#         -stats Processing/Wu61.merged.stats \
#         --ob-priors Processing/read-orientation-model.tar.gz \
#         -O Processing/Wu61.mutect2filtered.vcf