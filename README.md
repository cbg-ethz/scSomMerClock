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
