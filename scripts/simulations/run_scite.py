#!/usr/bin/env python3

import os
import re
import argparse
import subprocess
import numpy as np
from utils import get_sample_dict_from_vcf


def run_scite_subprocess(vcf_file, exe, steps, out_dir, prefix='', fd=0.001,
        ad=0.2, include='', exclude='', verbose=False):
    if not os.path.exists(exe):
        raise RuntimeError(
            'SCITE not compiled: run "g++ *.cpp -o SCITE -std=c++11" inside ' \
            f'{os.path.dirname(exe)}'
        )

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    data_raw, _ = get_sample_dict_from_vcf(vcf_file,
        GT=True, include=include, exclude=exclude)
    data_list = []
    for cell_data in data_raw.values():
        data_list.append([int(i) for i in cell_data])
    data = np.array(data_list).T

    if not prefix:
        try:
            run_no = re.search('\d{4}', os.path.basename(vcf_file)).group()
        except AttributeError:
            run_no = ''
        out_prefix = f'scite_tree.{run_no}'
    else:
        out_prefix = prefix

    data_file = os.path.join(out_dir, f'{out_prefix}.csv')
    np.savetxt(data_file, data.astype(int), delimiter=' ', newline='\n', fmt='%d')
    no_muts, no_cells = data.shape

    mut_file = os.path.join(out_dir, f'{prefix}muts.txt')
    if not os.path.exists(mut_file):
        with open(mut_file, 'w') as f_mut:
            f_mut.write('\n'.join([f'm{i}' for i in range(no_muts)]))

    out_files = os.path.join(out_dir, out_prefix)

    cmmd = ' '.join(
        [exe, '-i', data_file, '-transpose', '-r 1', '-n', str(no_muts),
        '-m', str(no_cells), '-l', str(steps), '-fd', str(fd), '-ad', str(ad),
        '-e 0.1', '-z', '-a',  '-o', out_files, '-names', mut_file,
        '-max_treelist_size 1']
    )

    if verbose:
        print(f'output directory:\n{out_dir}')
        print(f'\nShell command:\n{cmmd}\n')

    SCITE = subprocess.Popen(
        cmmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = SCITE.communicate()
    SCITE.wait()
    stdout = str(stdout)

    if stderr:
        for i in stdout.split('\\n'):
            print(i)
        raise RuntimeError(f'SCITE Error: {stderr}')

    with open(f'{out_files}.log', 'w') as f:
        f.write(stdout)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input vcf file')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory.')
    parser.add_argument('-p', '--prefix', type=str, default='',
        help='Output file prefixe.')
    parser.add_argument('-s', '--steps', type=int, default=500000,
        help='Number of mcmc steps. Default=500000.')
    parser.add_argument('-e', '--exe', type=str, help='Path to scite exe.')
    parser.add_argument('--FP', type=float, default=0.001,
        help='False positive rate. Default=0.001.')
    parser.add_argument('--FN', type=float, default=0.2,
        help='False negative rate. Default=0.2.')
    parser.add_argument('--verbose', action='store_true', help='Print to stdout')
    parser.add_argument('-ex', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test.')
    parser.add_argument('-in', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        out_dir, prefix_raw = os.path.split(snakemake.output[0])
        prefix = prefix_raw.replace('_ml0.newick', '')
        run_scite_subprocess(
            vcf_file=snakemake.input[0],
            exe=snakemake.params.exe,
            steps=snakemake.params.steps,
            out_dir=out_dir,
            prefix=prefix
        )
    else:
        args = parse_args()
        if not args.out_dir:
            args.out_dir = os.path.join(os.path.dirname(args.input), 'scite_dir')

        run_scite_subprocess(
            vcf_file=args.input,
            exe=args.exe,
            steps=args.steps,
            out_dir=args.out_dir,
            prefix=args.prefix,
            fd=args.FP,
            ad=args.FN,
            include=args.include,
            exclude=args.exclude,
            verbose=args.verbose
        )