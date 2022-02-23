#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import shutil
import subprocess


def run_bash(cmd):
    print(f'Running:\n{cmd}')
    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()

    print(str(stdout), str(stderr))


def run_plotting_from_folder(args):
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    for vcf_file_raw in os.listdir(args.in_dir):
        if not (vcf_file_raw.endswith('.vcf') or vcf_file_raw.endswith('.vcf.gz')):
            continue

        if vcf_file_raw.endswith('.vcf'):
            shutil.move(os.path.join(args.in_dir, vcf_file_raw),
                 os.path.join(args.in_dir, vcf_file_raw + '.gz'))
            vcf_file = os.path.join(args.in_dir, vcf_file_raw + '.gz')
        else:
            vcf_file = os.path.join(args.in_dir, vcf_file_raw)
        tree_file = vcf_file.replace('.vcf.gz', '.newick')

        assert os.path.exists(tree_file)

        file_name_raw = os.path.basename(vcf_file).replace('.vcf.gz', '')
        filters, tree = file_name_raw.split('_')[-2:]

        filter_dir = os.path.join(args.out_dir, filters)
        if not os.path.exists(filter_dir):
            os.makedirs(filter_dir)
        tree_dir = os.path.join(filter_dir, tree)
        if not os.path.exists(tree_dir):
            os.makedirs(tree_dir)

        out_file_raw = os.path.join(tree_dir, file_name_raw)
        out_file = out_file_raw + f'_w{args.w_max:.0f}_mapped.pdf'

        if os.path.exists(out_file) and not args.replace:
            print(f'Existing tree plot: {out_file}')
            continue

        cmd = f'python {args.exe} {vcf_file} {tree_file} -w {args.w_max} ' \
                f'-o {plot_file_raw} -b -p'
        run_bash(cmd)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_dir', type=str,  help='vcf master file')
    parser.add_argument('-o', '--out_dir', type=str, default = 'phylogeny_figs',
        help='Output directory.')
    parser.add_argument('-e', '--exe', type=str,
        default='simulations/scripts/get_poisson_tree_LRT.py',
        help='Poisson Tree exe.')
    parser.add_argument('-w', '--w_max', type=int, default=500,
        help='W_max value for Poisson Tree test.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    run_plotting_from_folder(args)