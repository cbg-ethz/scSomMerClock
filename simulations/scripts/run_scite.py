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
            'SCITE not compiled: run "g++ *.cpp -o SCITE -std=c++11" inside '
            '{}'.format(os.path.dirname(exe))
        )

    try:
        run_no = '.' + re.search('\d\d\d\d', os.path.basename(vcf_file)).group()
    except AttributeError:
        run_no = ''

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    data_raw, _ = get_sample_dict_from_vcf(vcf_file,
        GT=True, include=include, exclude=exclude)
    data_list = []
    for cell_data in data_raw.values():
        data_list.append([int(i) for i in cell_data])
    data = np.array(data_list).T

    data_file = os.path.join(out_dir, '{}SCITE{}.csv'.format(prefix, run_no))
    np.savetxt(data_file, data.astype(int), delimiter=' ', newline='\n', fmt='%d')
    no_muts, no_cells = data.shape

    mut_file = os.path.join(out_dir, '{}muts.txt'.format(prefix))
    if not os.path.exists(mut_file):
        with open(mut_file, 'w') as f_mut:
            f_mut.write('\n'.join(['m{}'.format(i) for i in range(no_muts)]))

    if not prefix:
        out_prefix = 'scite_tree{}'.format(run_no)
    else:
        out_prefix = prefix
    out_files = os.path.join(out_dir, out_prefix)

    cmmd = ' '.join(
        [exe, '-i', data_file, '-transpose', '-r 1', '-n', str(no_muts),
        '-m', str(no_cells), '-l', str(steps), '-fd', str(fd), '-ad', str(ad),
        '-e 0.1', '-z', '-a',  '-o', out_files, '-names', mut_file,
        '-max_treelist_size 1']
    )

    if verbose:
        print('output directory:\n{}'.format(out_dir))
        print('\nShell command:\n{}\n'.format(cmmd))

    SCITE = subprocess.Popen(
        cmmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = SCITE.communicate()
    SCITE.wait()
    stdout = str(stdout)

    if stderr:
        for i in stdout.split('\\n'):
            print(i)
        raise RuntimeError('SCITE Error: {}'.format(stderr))

    with open('{}.log'.format(out_files), 'w') as f:
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
        out_dir = os.path.sep.join(
            snakemake.input[0].split(os.path.sep)[:-2] + ['scite_dir'])
        run_scite_subprocess(
            vcf_file=snakemake.input[0],
            exe=snakemake.params.exe,
            steps=snakemake.params.steps,
            out_dir=out_dir,
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