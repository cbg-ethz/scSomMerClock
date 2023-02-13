# Single-Cell SOMatic MolEculaR Clock testing

This repository contains scripts for running
- the **Poisson Tree (PT) Test**
- and, in a subfolder (AnalysisPipelines), scripts for
  - the processing of real scDNA-seq data
  - the analysis of real scDNA-seq data
  - the simulation of scDNA-seq data (via coalescent)
  - the analysis and plotting of simulated scDNA-seq data

# Installation
## Requirements
- python3.X:
    - ete3
    - numpy
    - pandas
    - scipy

The requirements cant be installed using pip:
```bash
python -m pip install ete3 pandas scipy
```

# Usage
The **PT test** can be run with the following shell command:
```bash
python run_PT_test.py <VCF_FILE> <NEWICK_TREE_FILE> [-o] [-excl] [-incl] [-w] [-FN] [-FP]
```

## Input files
The **PT test** requires two input files:
- Called variants in VCF format ([VCF info](https://samtools.github.io/hts-specs/VCFv4.2.pdf)), where each sample is a cell
- An inferred phylogenetic tree in newick format (cell names need to be the same as in the VCF).

> ## Note
>  Trees can be inferred, for example, with [CellPhy](https://github.com/amkozlov/cellphy) or [infSCITE](https://github.com/cbg-ethz/infSCITE); both outputs are compatible with the **PT test**

## Optional Arguments
- `-o <str>`, Output file. Default = <VCF_FILE>.poissonTree_LRT.tsv.
- `-excl <str>`, Regex pattern for samples/cells to exclude. Default = none.
- `-incl <str>`, Regex pattern for samples/cells to include. If set, only these samples/cells are included. Default = all cells.
- `-w <list of int>`, Maximum weight values. Defaut = 100, 200, ..., 1000'.
- `-FN <float>`, Estimated FN rate (for CellPhy and infSCITE: inferred from .log/stdout file).
- `-FP <float>`, Estimated FP rate (for CellPhy and infSCITE: inferred from .log/stdout file).

# Example
To run the PT test on the simulated data in the  `example_data` folder, execute
```bash
python run_PT_test.py example_data example_data/data_simulated_clock.vcf.gz example_data/data_simulated_clock.raxml.bestTree
```
or
```bash
python run_PT_test.py example_data example_data/data_simulated_noclock.vcf.gz example_data/data_simulated_noclock.raxml.bestTree
```

The former data is simulated under a molecular clock, the later with a deviation from the clock (evolutionary rate amplified by 5x in a subtree)
> ## Note
>  FN and FP rate are inferred from the `.raxml.log` file
