#!/usr/bin/env python3

import os
import re
from statistics import mean, stdev
import numpy as np
from scipy.stats.distributions import chi2


def get_muts_per_cell(vcf_file, exclude):
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
    for i in exclude_i:
        samples.pop(i)

    return samples


def test_poisson(in_files, out_file, exclude='', alpha=0.05):
    avg = [[], [], 0, 0]
    out_str = ''
    for in_file in in_files:
        run = os.path.basename(in_file).split('.')[1]
        muts = np.array(get_muts_per_cell(in_file, exclude))

        mean_muts = muts.mean()

        LR = 2 * np.sum(muts * np.log(muts / mean_muts))
        # LR2 =  np.sum((muts - mean_muts)**2) / mean_muts
        avg[0].append(LR)
        p_val = chi2.sf(LR, len(muts) - 1)
        avg[1].append(p_val)
        if p_val < alpha:
            hyp = 'H1'
            avg[2] += 1
        else:
            hyp = 'H0'
            avg[3] += 1

        out_str += f'{run}\t{LR:0>5.2f}\t{p_val:.2E}\t{hyp}\n'

    avg_line = f'\nAvg.\t{mean(avg[0]):0>5.2f}\t{mean(avg[1]):.2E}\t' \
            f'H0:{avg[3]};H1:{avg[2]}\n'
    with open(out_file, 'w') as f_out:
        f_out.write('run\t-2logLR\tp-value\thypothesis\n')
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        test_poisson(snakemake.input, snakemake.output[0],
            snakemake.params.exclude)
    else:
        import argparse
        args = parse_args()
        test_poisson(args.input, args.output, args.exclude, args.alpha)