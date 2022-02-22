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
    w_max = 500

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

        print(vcf_file)
        assert os.path.exists(tree_file)

        plot_file_raw = vcf_file.replace('.vcf', '')
        plot_file = plot_file_raw + f'_w{w_max:.0f}_mapped.png'

        cmd = f'python {args.exe} {vcf_file} {tree_file} -w {w_max} ' \
                f'-o {plot_file_raw} -b -p'
        run_bash(cmd)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_dir', type=str,  help='vcf master file')
    parser.add_argument('-o', '--out_dir', type=str, default = '',
        help='Output directory.')
    parser.add_argument('-e', '--exe', type=str,
        default='simulations/scripts/get_poisson_tree_LRT.py',
        help='Poisson Tree exe.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    run_plotting_from_folder(args)