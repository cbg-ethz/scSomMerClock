# MolClockAnalysis
ScDNAseq data processing pipeline for molecular clock analysis


## Requirements
- cutadapt >1.18
- gatk >4.0
- gcc >6.4
- gcccore >6.4.0
- numpy (python 2.7)
- picard >2.18
- python >3.7
- python 2.7
- pysam
- samtools
- snakemake >5.4.5


## Folder Structure
Data directory:
```bash
<DATA_DIR>
 \_ <CELLNAMES>.txt
 \_ 
 \_ ...
```
Resources directory:
```bash
<RESOURCES_DIR>
 \_ <REF>.fa
 \_ <REF>.dict
 \_ <DBSNP>.vcf
 \_ <DBSNP>.vcf.idx
 \_ <INDEL1>.vcf
 \_ <INDEL1>.vcf.idx
 \_ <INDEL2>.vcf
 \_ <INDEL2>.vcf.idx
 \_ <EXON_TARGET>.bed
 \_ ...
```

> **Note**: To create .idx files, use: ``gatk IndexFeatureFile -F <FILE>.vcf``

## How to run
Load snakemake (e.g. via module or conda):
```bash
module load snakemake
# conda activate <SNAKEMAKE_ENV>
```

Run on a hpc cluster from the base directory:
```bash
bash ./wrapper.sh -c configs/<DATASET> --profile hpc/<slurm|lsf>
```



## How to run
Enter the following to the command line:
```bash
snakemake -j 1 -s Snakefile_sim --configfile config.simulations_clock.yaml --restart-times=0
```

To run on CESGA hpc cluster
```bash
module load snakemake
snakemake -j 198 -s Snakefile_sim --configfile configs/config.simulations_clock_scWGA_NGS.yaml --use-envmodules -k --profile ../hpc/slurm --scheduler greedy &> logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out &
```


To run on ETHZ hpc cluster, you need a conda environment with snakemake installed. If it is called snakemake, run:
```bash
conda activate snakemake
snakemake -s Snakefile_sim --configfile configs/config.simulations_clock_scWGA_NGS.yaml --use-conda -k --profile ../hpc/lsf &> logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out &
```