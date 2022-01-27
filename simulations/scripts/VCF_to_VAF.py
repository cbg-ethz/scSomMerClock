#!/usr/bin/env python3

import argparse
import gzip
import os
import re

import numpy as np


MUT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
MUT_REV = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
MUT_VAR = {'A': np.array([1, 2, 3]), 'C': np.array([0, 2, 3]),
    'G': np.array([0, 1, 3]), 'T': np.array([0, 1, 2])}

def vcf_to_vaf(vcf_file, incl_re='', excl_re='healthycell', fmin=1e-6, fmax=1):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude = []
    vafs = []
    # Compatible with mobster dN/dS calculation
    out_str = 'chrom\tfrom\tref\talt\tDP\talt_count\tVAF'

    with file_stream as f_in:
        for line in f_in:
            try:
                line = line.decode()
            except:
                pass
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    if excl_re != '':
                        print('Exluded:')
                        for s_i, sample in enumerate(sample_names):
                            if re.fullmatch(excl_re, sample):
                                exclude.append(s_i)
                                print('\t{}'.format(sample))
                        if len(exclude) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(excl_re))
                    if incl_re != '':
                        print('Include:')
                        for s_i, sample in enumerate(sample_names):
                            if not re.fullmatch(incl_re, sample):
                                exclude.append(s_i)
                            else:
                                if not re.fullmatch(excl_re, sample):
                                    print('\t{}'.format(sample))
                        if len(exclude) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(incl_re))
                continue
            elif line.strip() == '':
                continue

            line_cols = line.strip().split('\t')
            elem_ids = line_cols[8].split(':')

            reads = np.zeros(4, dtype=float) # Order: ACGT

            for s_i, s_rec in enumerate(line_cols[9:]):
                if s_i in exclude:
                    continue
                s_rec_elms = s_rec.split(':')
                try:
                    s_read = s_rec_elms[read_id].split(',')
                except NameError:
                    read_id = elem_ids.index('RC')
                    s_read = s_rec_elms[read_id].split(',')

                reads += np.array(s_read, dtype=float)

            var_all = MUT_VAR[line_cols[3]]
            DP = reads.sum()

            for alt_id, alt_count in enumerate(reads):
                if alt_id == MUT[line_cols[3]]:
                    continue
                alt = MUT_REV[alt_id]
                vaf = alt_count / DP
                # chr pos DP alt_count VAF
                if vaf >= fmin and vaf < fmax:
                    out_str += f'\n{line_cols[0]}\t{line_cols[1]}\t' \
                        f'{line_cols[3]}\t{alt}\t{DP}\t{alt_count}\t{vaf:.6f}'
                    vafs.append(vaf)

    return np.array(vafs), out_str


def show_VAF_dist(arr):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.hist(arr, bins=100, range= (0, 1))
    ax.set_ylabel(f'counts (n={arr.size})', fontsize=20)
    ax.set_xlabel(f'VAF', fontsize=20)
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92, hspace=0.5)
    plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in VCF format.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file. Default: <VCF>.VAF')
    parser.add_argument('-min', '--fmin', type=float, default=1e-6,
        help='Minimum VAF value: <= are excluded. Default: 0.001')
    parser.add_argument('-max', '--fmax', type=float, default=1.0,
        help='Maximum VAF value: > are excluded. Default: 1.0')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include.')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        vafs, out_str = vcf_to_vaf(snakemake.input.vcf)
        out_file = snakemake.output.vaf
    else:
        args = parse_args()
        vafs, out_str = vcf_to_vaf(args.vcf, args.include, args.exclude, args.fmin, args.fmax)

        if args.output:
            out_file = args.output
        else:
            out_file = args.vcf + '.VAF'

    with open(out_file, 'w') as f:
        f.write(out_str)

    if 'snakemake' in globals():
        os.remove(snakemake.input.vcf)
    # show_VAF_dist(vafs)