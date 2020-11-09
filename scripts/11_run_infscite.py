#!/usr/bin/env python3

import argparse
import os
import subprocess
from tempfile import NamedTemporaryFile


def parse_args():
    parser = argparse.ArgumentParser(
        description='*** Run inf-sSCITE on genotype matrix ***')
    parser.add_argument('input', type=str,
        help='Genotype matrix where columns are samples and rows are mutations.')
    parser.add_argument('-e', '--exe', type=str, required=True,
        help='Path to inf-scite exe.')
    parser.add_argument('-o', '--outdir', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    parser.add_argument('-r', '--repetitions', type=int, default=1,
        help='Desired number of repetitions of the MCMC. Default = 1.')
    parser.add_argument('-l', '--length', type=int, default=100000,
        help='MCMC chain length. Default = 100000.')
    parser.add_argument('-er', '--error', type=float, default=0.2,
        help='Error learning rate. Default = 0.2.')
    parser.add_argument('-ad1', type=float, default=0.25,
        help='Estimated rate of REF/ALT called as REF/REF. Default = 0.25.')
    parser.add_argument('-ad2', type=float, default=0.25,
        help='Estimated rate of REF/ALT called as ALT/ALT. Default = 0.25.')
    parser.add_argument('-fd', type=float, default=0.0001,
        help='Estimated rate of REF/REF called as REF/ALT. Default = 0.0001.')
    parser.add_argument('-cc', type=float, default=0.000001,
        help='Estimated rate of REf/REF called as ALT/ALT. Default = 0.000001.')
    parser.add_argument('-pf', '--permanent_files', action='store_true',
        help='Write input files as permanent, not temp files.')
    args = parser.parse_args()
    return args


def load_data(in_file):
    # Read in all lines from in_file
    with open(in_file, 'r') as f:
        lines = f.read().strip().split('\n')

    # Get separator
    seps = ['\t', ' ', ',']
    sep_counts = [0, 0, 0]
    for i, sep in enumerate(seps):
        sep_counts[i] = lines[0].count(sep)
    sep = seps[sep_counts.index(max(sep_counts))]

    header_row = False
    header_line = ''
    for el in lines[0].split(sep):
        try:
            el_float = float(el)
        except ValueError:
            if el == ' ':
                continue
            header_row = True
            header_line = lines.pop(0)
            break    
        else:
            if el_float not in [0, 1, 2, 3]:
                header_row = True
                header_line = lines.pop(0)
                break
    
    index_col = False
    # Check first 5 lines if index column is set
    for i, line in enumerate(lines[:5]):
        first_el = line.split(sep)[0]
        try:
            first_el_flt = float(first_el)
        except ValueError:
            if first_el == ' ':
                continue
            index_col = True
            break
        else:
            if first_el_flt not in [0, 1, 2, 3]:
                index_col = True
                break

    samples = [i for i in header_line.split(sep) if i]
    if index_col:
        samples.pop(0)
    index = []
    data = []
    for line in lines:
        line_el = line.split(sep)
        if index_col:
            index.append(line_el.pop(0))
        data.append(line_el)
        
    return data, samples, index


def get_data_infSCITE(args):
    data, samples, index = load_data(args.input)
    scite_args = {'n': len(data), 'm': len(data[0]), 'names': '', 'samples': ''}

    # Write data to file
    data_out = [' '.join(i) for i in data]
    if args.permanent_files:
        data_file = os.path.join(args.outdir, 'infSCITE_genotype_matrix.csv')
        with open(data_file, 'w') as data_f:
            data_f.write('\n'.join(data_out))
    else:
        with NamedTemporaryFile(delete=False, mode='w') as data_f:
            data_f.write('\n'.join(data_out))
        data_file = data_f.name
    scite_args['input'] = data_file
    # Write mutation names to file if given
    if index:
        if args.permanent_files:
            mut_file = os.path.join(args.outdir, 'infSCITE_mutations.txt')
            with open(mut_file, 'w') as mut_f:
                mut_f.write('\n'.join(index))
        else:
            with NamedTemporaryFile(delete=False, mode='w') as mut_f:
                mut_f.write('\n'.join(index))
            mut_file = mut_f.name
        scite_args['names'] = '-names {}'.format(mut_file)
    # Write sample names to file if given
    if samples:
        if args.permanent_files:
            sample_file = os.path.join(args.outdir, 'infSCITE_samples.txt')
            with open(sample_file, 'w') as sample_f:
                sample_f.write('\n'.join(samples))
        else:
            with NamedTemporaryFile(delete=False, mode='w') as sample_f:
                sample_f.write('\n'.join(samples))
            sample_file = sample_f.name
        scite_args['samples'] = '-samples {}'.format(sample_file)
    
    return scite_args


def main(args):
    if args.outdir == '':
        args.outdir = os.path.dirname(args.input)

    scite_args = get_data_infSCITE(args)

    os.makedirs(args.outdir, exist_ok=True)
    out_tree = os.path.join(args.outdir, 'infSCITE_tree')

    run_scite = '{exe} -i {input} -r {r} -l {l} -n {n} -m {m} -fd {fd} ' \
            '-ad {ad1} {ad2} -cc {cc} -e {e} -transpose -o {o} {names} {samples}' \
        .format(exe=args.exe, r=args.repetitions, l=args.length, 
            fd=args.fd, ad1=args.ad1, ad2=args.ad2, cc=args.cc, e=args.error,
            o=out_tree, **scite_args)
    print('\nRunning:\n{}\n'.format(run_scite))

    process = subprocess.Popen(run_scite,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    stdout = process.stdout.read().decode('utf-8').strip()

    print('Stdout:{}\n'.format(stdout))


if __name__ == '__main__':
    args = parse_args()
    main(args)