#!/usr/bin/env python3

import os
import re
from statistics import mean, stdev
import numpy as np
from scipy.stats.distributions import chi2


def get_muts_per_cell(vcf_file, exclude, include):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude_i = []
    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    samples = [0 for i in range(len(sample_names))]
                    if exclude != '':
                        print('Exluded:')
                        for s_i, sample in enumerate(sample_names):
                            if re.fullmatch(exclude, sample):
                                exclude_i.append(s_i)
                                print('\t{}'.format(sample))
                        if len(exclude_i) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(exclude))
                    if include != '':
                        print('Include:')
                        for s_i, sample in enumerate(sample_names):
                            if not re.fullmatch(include, sample):
                                exclude_i.append(s_i)
                            else:
                                if not re.fullmatch(exclude, sample):
                                    print('\t{}'.format(sample))
                        if len(exclude_i) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(include))
                continue

            elif line.strip() == '':
                continue

            # VCF records
            line_cols = line.strip().split('\t')
            bases = [line_cols[3]] + line_cols[4].split(',')
            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
                    continue

                s_rec_ref = gt[0]
                s_rec_alt = gt[-1]
                if s_rec_ref == '.' or s_rec_alt == '.':
                    continue

                if s_rec_alt != '0' or s_rec_ref != '0':
                    samples[s_i] += 1
                    
                # gt_pred = '{}|{}' \
                #     .format(bases[int(s_rec_ref)], bases[int(s_rec_alt)])
                # gt_true = s_rec[-3:]
                # if gt_pred != gt_true and gt_pred != gt_true[::-1]:
                #     print('True GT: {}; inferred GT: {}'.format(gt_true, gt_pred))

    # Remove excluded samples
    for i in sorted(exclude_i, reverse=True):
        samples.pop(i)
    # Sanity check: remove sample without name
    try:
        samples.pop(sample_names.index(''))
    except (ValueError, KeyError):
        pass
    
    return samples


def test_poisson(in_files, out_file, exclude='', include='', alpha=0.05):
    avg = [[], [], [], [], 0, 0]
    out_str = ''
    for in_file in in_files:
        run = os.path.basename(in_file).split('.')[1]
        muts = np.array(get_muts_per_cell(in_file, exclude, include))

        mean_muts = muts.mean()

        ll_clock = -np.sum(muts * np.log(mean_muts) - mean_muts)
        ll_uncon = -np.sum(muts * np.log(muts) - muts)
        LR = 2 * (ll_clock - ll_uncon)
        # LR2 =  np.sum((muts - mean_muts)**2) / mean_muts
        dof = len(muts) - 1
        p_val = chi2.sf(LR, dof)

        for i, j in [(0, ll_clock), (1, ll_uncon), (2, LR), (3, p_val),]:
            avg[i].append(j)

        if p_val < alpha:
            hyp = 'H1'
            avg[5] += 1
        else:
            hyp = 'H0'
            avg[4] += 1

        out_str += f'{run}\t{ll_clock:0>5.2f}\t{ll_uncon:0>5.2f}\t{LR:0>5.2f}\t' \
            f'{dof}\t{p_val:.2E}\t{hyp}\n'

    avg_line = f'\nAvg.\t{mean(avg[0]):0>5.2f}\t{mean(avg[1]):0>5.2f}\t' \
            f'{mean(avg[2]):0>5.2f}\t\t{mean(avg[3]):0>2.5f}\tH0:{avg[4]};H1:{avg[5]}\n'
    with open(out_file, 'w') as f_out:
        f_out.write('run\tH0\tH1\t-2logLR\tdof\tp-value\thypothesis\n')
        f_out.write(out_str + avg_line)

    print(avg_line)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
        help='Significance threshold. Default = 0.05.')
    parser.add_argument('-e', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-i', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''
        test_poisson(snakemake.input, snakemake.output[0],
            snakemake.params.exclude, snakemake.params.include)
    else:
        import argparse
        args = parse_args()
        test_poisson(args.input, args.output, args.exclude, args.include,
            args.alpha)