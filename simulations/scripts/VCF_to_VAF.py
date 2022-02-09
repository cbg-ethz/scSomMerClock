#!/usr/bin/env python3

import argparse
import gzip
import os
import re

import numpy as np


MUT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
MUT_REV = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
# MUT_VAR = {'A': np.array([1, 2, 3]), 'C': np.array([0, 2, 3]),
#     'G': np.array([0, 1, 3]), 'T': np.array([0, 1, 2])}

def vcf_to_vaf_sc(vcf_file, incl_re='', excl_re='healthycell', fmin=1e-6, fmax=1):

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

            for alt_id in np.argsort(reads)[::-1][:2]:
                if alt_id == MUT[line_cols[3]]:
                    continue
                alt_count = reads[alt_id]
                if alt_count <= 2:
                    continue

                DP = reads.sum()
                alt_count = reads[alt_id]
                alt = MUT_REV[alt_id]
                vaf = alt_count / DP

                if vaf >= fmin and vaf < fmax:
                    out_str += f'\n{line_cols[0]}\t{line_cols[1]}\t' \
                        f'{line_cols[3]}\t{alt}\t{DP}\t{alt_count}\t{vaf:.6f}'
                    vafs.append(vaf)

    # show_VAF_dist(vafs)
    return np.array(vafs), out_str


def vcf_to_vaf_bulk(vcf_file, fmin=1e-6, fmax=1, dmin=10):

    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

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
                    vafs = [[] for i in sample_names]
                    out_str = ['chrom\tfrom\tref\talt\tDP\talt_count\tVAF'] \
                        * len(sample_names)
                continue
            elif line.strip() == '':
                continue

            line_cols = line.strip().split('\t')
            elem_ids = line_cols[8].split(':')

            # Skip indels
            if len(line_cols[3]) > 1 or len(line_cols[4]) > 1:
                continue

            for s_i, s_rec in enumerate(line_cols[9:]):
                s_rec_elms = s_rec.split(':')

                try:
                    alt_count = s_rec_elms[ad_id].split(',')[1]
                    DP = int(s_rec_elms[depth_id])
                    vaf = float(s_rec_elms[vaf_id])
                except NameError:
                    ad_id = elem_ids.index('AD')
                    depth_id = elem_ids.index('DP')
                    vaf_id = elem_ids.index('AF')

                    alt_count = s_rec_elms[ad_id].split(',')[1]
                    DP = int(s_rec_elms[depth_id])
                    vaf = float(s_rec_elms[vaf_id])

                # Skip wt
                if s_rec_elms[0] == '0/0' or s_rec_elms[0] == '0|0':
                    continue

                if vaf >= fmin and vaf < fmax and DP >= dmin:
                    out_str[s_i] += f'\n{line_cols[0]}\t{line_cols[1]}\t' \
                        f'{line_cols[3]}\t{line_cols[4]}\t{DP}\t{alt_count}\t' \
                        f'{vaf:.6f}'
                    vafs[s_i].append(vaf)

    excl = [i for i, j in enumerate(vafs) if len(j) == 0]
    vafs = [j for i, j in enumerate(vafs) if i not in excl]
    out_str = [j for i, j in enumerate(out_str) if i not in excl]

    return vafs, zip(sample_names, out_str)


def show_VAF_dist(arr):
    import matplotlib.pyplot as plt

    if isinstance(arr, list):
        arr = np.array(arr)

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.hist(arr, bins=100, range=(0, 1))
    ax.set_ylabel(f'counts (n={arr.size})', fontsize=20)
    ax.set_xlabel(f'VAF', fontsize=20)
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92, hspace=0.5)
    plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in VCF format.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file/dir. Default: <VCF>.VAF')
    parser.add_argument('-min', '--fmin', type=float, default=1e-6,
        help='Minimum VAF value: <= are excluded. Default: 0.001')
    parser.add_argument('-max', '--fmax', type=float, default=1.0,
        help='Maximum VAF value: > are excluded. Default: 1.0')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include.')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude.')
    parser.add_argument('-b', '--bulk', action='store_true',
        help='If set, assumes mutect2 VCF output.')
    parser.add_argument('-dmin', '--depthmin', type=int, default=10,
        help='Min. depth for bulk/biological data.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        vafs, out_str = vcf_to_vaf_sc(snakemake.input.vcf)
        with open(snakemake.output.vaf, 'w') as f:
            f.write(out_str)
    else:
        args = parse_args()
        if args.output:
            out_file = args.output
        else:
            out_file = args.vcf

        if args.bulk:
            vafs, out_data = vcf_to_vaf_bulk(args.vcf, args.fmin, args.fmax,
                args.depthmin)

            for sample_name, out_str in out_data:
                out_file_s = f'{out_file}_{sample_name}.VAF'
                with open(out_file_s, 'w') as f:
                    f.write(out_str)
        else:
            vafs, out_str = vcf_to_vaf_sc(args.vcf, args.include, args.exclude,
                args.fmin, args.fmax)

            with open(out_file, 'w') as f:
                f.write(out_str)

    # show_VAF_dist(vafs)