## How to run
Enter the following to the command line:
```bash
snakemake -j 1 -s Snakefile_sim --configfile config.simulations_clock.yaml
```

To run on CESGA hpc cluster
```bash
snakemake -j 198 -s Snakefile_sim --configfile config.simulations_clock.yaml -k --profile ../hpc/slurm &> logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out
```