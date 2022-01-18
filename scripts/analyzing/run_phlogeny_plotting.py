#!/usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
import subprocess

base_dir = '/home/uvi/be/nbo/data/data/'

data_dirs = {
    'H65_Monica': ['all', 'cancer', 'normal'],
    'Li55': ['all', 'cancer', 'normal'],
    'Lo-P1': ['all'],
    'Lo-P2': ['all'],
    'Lo-P3': ['all'],
    'Ni8_Monica': ['all', 'cancer'],
    'S21_P1': ['all'],
    'S21_P2': ['all', 'cancer', 'left'],
    'W32_Monica': ['all', 'aneuploid', 'cancer', 'haploid', 'normal'],
    'W55': ['all', 'cancer', 'normal'],
    'Wu61': ['all', 'cancer', 'cancer_C', 'cancer_CA', 'cancer_C_CA', 'normal',
        'polyps'],
    'Wu63': ['all', 'cancer', 'cancer_polyps', 'normal', 'polyps'],
    'X25': ['all', 'cancer', 'normal']
}
data_filters = ['all', '33nanFilter', '50nanFilter']

MODULE_STR = 'module load ete;'


def run_bash(cmd_raw, bsub=True):
    if bsub:
        cmd = f"sbatch -t 30 -p amd-shared --qos amd-shared --mem 2G " \
            f"--wrap '{MODULE_STR} {cmd_raw}'"
    else:
        cmd = cmd_raw

    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()

    print(f'Running: {cmd}')
    if not bsub:
        print(str(stdout), str(stderr))
    print('\n')


def run_phlogeny_plotting(args):
    for data_dir, sub_dirs in data_dirs.items():
        data_set = data_dir.replace('_Monica', '')
        # Iterate sub sets
        for sub_dir in sub_dirs:
            vcf_dir = os.path.join(base_dir, data_dir, 'ClockTest', sub_dir)
            # Iterate nan filters
            for filters in data_filters:
                vcf_name = f'{data_set}.{filters}.vcf.gz'
                vcf_file = os.path.join(vcf_dir, vcf_name)

                cellphy = ('cellphy', vcf_file + '.raxml.bestTreemapped.newick')
                scite = ('scite', os.path.join(vcf_dir, 'scite_dir',
                    f'{data_set}.{filters}_ml0.newickmapped.newick'))

                for tree, in_file in [scite, cellphy]:
                    if not os.path.exists(in_file):
                        print(f'Missing Tree file: {in_file}')
                        continue

                    out_file = os.path.join(args.out_dir,
                        f'Phylogeny_{tree}_{data_set}_{sub_dir}_{filters}.pdf')
                    if not os.path.exists(out_file) or args.replace:
                        cmd =  f'python {args.exe} -i {in_file} -o {out_file}'
                        run_bash(cmd, args.local)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--out_dir', type=str,
        default='poisson_tests_all/phylogenies', help='Output directory.')
    parser.add_argument('-e', '--exe', type=str,
        default='simulations/plotting/plot_tree.py', help='Tree plotting script.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    run_phlogeny_plotting(args)




