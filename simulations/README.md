## How to run
Enter the following to the command line:
```bash
snakemake -j 1 -s Snakefile_sim --configfile config.simulations_clock.yaml
```

To run on CESGA hpc cluster
```bash
module load snakemake
snakemake -j 198 -s Snakefile_sim --configfile config.simulations_clock.yaml -k --profile ../hpc/slurm &> logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out
```


To run on ETHZ hpc cluster
```bash
source activate crp
snakemake -s Snakefile_sim --configfile config.simulations_clock.yaml --use-conda -k --profile ../hpc/lsf &> logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out
```