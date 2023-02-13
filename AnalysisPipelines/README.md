# Single-Cell SOMatic MolEculaR Clock testing - Analysis Pipelines

This folder contains the software for:
- the processing of real scDNA-seq data
- the analysis of real scDNA-seq data
- the simulation of scDNA-seq (via coalscent)
- the analysis and plotting of simulated scDNA-seq data

## Requirements

### general
- snakemake >5.4.5

### processing
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

### simulations, analysis & plotting

- python3.X:
    - ete3
    - matplotlib
    - numpy
    - pandas
    - scipy
    - seaborn
- R:
    - argparser
    - ggplot2
    - mobster
    - neutralitytestr
- cellcoal
- PAUP*
- CellPhy
- infSCITE/SCITE

## How to run scDNA-seq processing
Load snakemake (e.g. via module or conda):
```bash
module load snakemake
# conda activate <SNAKEMAKE_ENV>
```

Execute snakemake from the base directory:
```bash
snakemake -j 1 -s Snakefile --configfile config/data/<DATASET>.yaml
```


## How to run simulations
Enter the following to the command line:
```bash
snakemake -j 1 -s Snakefile_sim --configfile config/simulations/<SIM_SETUP>.yaml
```

To run on a hpc cluster with either slurm or lsf:
```bash
snakemake -j 200 -s Snakefile_sim --configfile configs/simulations/<SIM_SETUP>.yaml --use-conda -k --profile ../hpc/<lsf|slurm> &> logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out &
```

> **Note**: replace ``--use-conda`` with ``--use-envmodules`` if working with module