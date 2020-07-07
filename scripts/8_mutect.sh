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
                            bams_in+="-I $1 " 
                            normal_in="-normal $1 "
                            ;;
        -t | --tumor )      shift
                            bams_in+="-I $1 " 
                            ;;
    esac
    shift
done

[[ -z "$chr" ]] && { echo "Error: Chromosome not set"; exit 1; }
[[ -z "$REF" ]] && { echo "Error: Reference not set"; exit 1; }
[[ -z "$GERM_RES" ]] && { echo "Error: Germline resources not set"; exit 1; }
[[ -z "$PON" ]] && { echo "Error: Panel-Of-Normals not set"; exit 1; }

tumors_in=$(awk -v chr=$SLURM_ARRAY_TASK_ID '{print "-I "Processing"/"$0".real."chr".bam"}' bulk-tumor-name.txt | tr '\n' ' ')
bulk_in=$(awk -v chr=$SLURM_ARRAY_TASK_ID '{print "-I "Processing"/"$0".real."chr".bam"}' bulk-normal-name.txt | tr '\n' ' ')
bulk_name=$(cat bulk-normal-name.txt)

gatk Mutect2 \
    -R ${REF} \
    ${bams_in} \
    ${normal_in} \
    --genotype-pon-sites true \
    --germline-resource ${GERM_RES} \
    --panel-of-normals ${PON} \
    -O Calls/${chr}.mutect.vcf \
    --f1r2-tar-gz Calls/${chr}.f1r2.mutect.tar.gz
