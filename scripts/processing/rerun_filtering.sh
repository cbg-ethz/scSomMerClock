#!/bin/sh

module load pysam
scriptDir='/home/uvi/be/nbo/MolClockAnalysis/scripts'
DP=10
GQ=13
runCellphy=false

while [ "$1" != "" ]; do
    key=$1
    case ${key} in
        -i | --in)          shift
                            dataDir=$1
                            ;;
        -b | --bulk )       shift
                            bulk=$1
                            ;;
        -d | --DP)          shift
                            DP=$1
                            ;;
        -q | --GQ )         shift
                            GQ=$1
                            ;;
        -s | --scripts)     shift
                            scriptDir=$1
                            ;;
        -c | --cellphy)     runCellphy=true
                            ;;
    esac
    shift
done

set -Eeuxo pipefail
cores=$(nproc)

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 'X' 'Y'; do
    python ${scriptDir}/processing/10_summarize_vcf.py ${dataDir}/all.${chr}.vcf.gz -o ${dataDir} -bt ${bulk} -q ${GQ} -r ${DP}
done;

input=$(ls ${dataDir}/all_filtered.*.vcf);
python ${scriptDir}/processing/11_merge_filetered_vcfs.py $input -o ${dataDir};
if [ "$runCellphy" = true ]; then
    sbatch -t 1000 -p amd-shared --qos amd-shared --mem 4G --wrap "/home/uvi/be/nbo/cellphy/cellphy.sh -r -t ${cores} -z -l ${dataDir}/all_filtered.vcf";
fi
