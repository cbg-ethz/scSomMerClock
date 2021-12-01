#!/usr/bin/env python3

import gzip
import numpy as np
import os
import re
from statistics import mean, stdev
from scipy.stats.distributions import chi2
from tqdm import tqdm

LAMBDA_MIN = 1e-12


def get_muts_per_cell(vcf_file, exclude, include):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude_i = []
    with file_stream as f_in:
        for line in f_in:
            try:
                line = line.decode()
            except AttributeError:
                pass
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    sample_no = len(sample_names)
                    samples = np.zeros(sample_no)
                    include_i = np.full(sample_no, True)

                    if exclude != '':
                        print('Exluded:')
                        for s_i, sample in enumerate(sample_names):
                            if re.fullmatch(exclude, sample):
                                include_i[s_i] = False
                                print('\t{}'.format(sample))
                        if include_i.sum() == sample_no:
                            print('\nWARNING: no samples with pattern ' \
                                f'{exclude} in vcf!\n')
                    if include != '':
                        print('Include:')
                        for s_i, sample in enumerate(sample_names):
                            if not re.fullmatch(include, sample):
                                include_i[s_i] = False
                            else:
                                if not re.fullmatch(exclude, sample):
                                    print('\t{}'.format(sample))
                        if include_i.sum() == sample_no:
                            print('\nWARNING: no samples with pattern ' \
                                f' {include} in vcf!\n')
                    # Sanity check: remove sample without name
                    for s_i, sample in enumerate(sample_names):
                        if sample == '':
                            include_i[s_i] = False
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
    return samples[include_i]


def test_poisson_simulation(in_files, out_file,
            exclude='', include='', alpha=0.05):
    clock = re.search('clock(\d+.?\d*)_', out_file).group(1) == '0'
    exclude += '|healthycell'

    avg = [[], [], [], [], 0]
    out_str = ''
    for file_no, in_file in enumerate(sorted(in_files)):
        run = int(os.path.basename(in_file).split('.')[1])
        muts = np.array(get_muts_per_cell(in_file, exclude, include))
        mean_muts = muts.mean()

        h0 = -np.nansum(muts * np.log(mean_muts) - mean_muts)
        h1 = -np.nansum(muts * np.log(np.clip(muts, LAMBDA_MIN, None)) - muts)
        LR = 2 * (h0 - h1)
        dof = muts.size - 1
        p_val = chi2.sf(LR, dof)

        for i, j in [(0, h0), (1, h1), (2, LR), (3, p_val),]:
            avg[i].append(j)

        if p_val < alpha:
            hyp = 'H1'
            if not clock:
                avg[4] += 1
        else:
            hyp = 'H0'
            if clock:
                avg[4] += 1

        out_str += f'{run}\t{h0:0>5.2f}\t{h1:0>5.2f}\t{LR:0>5.2f}\t{dof}\t' \
            f'{p_val:.2E}\t{hyp}\n'

    avg_line = f'\n-1\t{mean(avg[0]):0>5.2f}\t{mean(avg[1]):0>5.2f}\t' \
            f'{mean(avg[2]):0>5.2f}\t{dof}\t{mean(avg[3]):0>2.5f}\t' \
            f'{avg[4]}/{file_no+1}\n'
    with open(out_file, 'w') as f_out:
        f_out.write('run\tH0\tH1\t-2logLR\tdof\tp-value\thypothesis\n')
        f_out.write(out_str.rstrip() + avg_line)

    print(avg_line)


def test_poisson_biological(in_files, out_file, excl='', incl='', alpha=0.05):
    out_str = ''
    for in_file in tqdm(sorted(in_files)):
        path_strs = in_file.split(os.path.sep)
        clock_dir_no = [i for i, j in enumerate(path_strs) if j == 'ClockTest'][0]
        dataset = path_strs[clock_dir_no - 1]
        subset = path_strs[clock_dir_no + 1]

        muts = np.array(get_muts_per_cell(in_file, excl, incl))
        mean_muts = muts.mean()

        h0 = -np.nansum(muts * np.log(mean_muts) - mean_muts)
        h1 = -np.nansum(muts * np.log(np.clip(muts, LAMBDA_MIN, None)) - muts)
        LR = 2 * (h0 - h1)
        dof = muts.size - 1
        p_val = chi2.sf(LR, dof)

        if p_val < alpha:
            hyp = 'H1'
        else:
            hyp = 'H0'

        out_str += f'{dataset}\t{subset}\t{h0:0>5.2f}\t{h1:0>5.2f}\t{LR:0>5.2f}' \
            f'\t{dof}\t{p_val:.2E}\t{hyp}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('dataset\tsubset\tH0\tH1\t-2logLR\tdof\tp-value\thypothesis\n')
        f_out.write(out_str.rstrip())


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
    parser.add_argument('-b', '--biological_data', action='store_true',
        help='Test true data (isntead of simulation data).')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''

        test_poisson_simulation(
            in_files=snakemake.input,
            out_file=snakemake.output[0],
            exclude=snakemake.params.exclude,
            include=snakemake.params.include,
        )
    else:
        import argparse
        args = parse_args()
        if args.biological_data:
            test_poisson_biological(
                in_files=args.input,
                out_file=args.output,
                excl=args.exclude,
                incl=args.include,
                alpha=args.alpha
            )
        else:
            test_poisson_simulation(
                in_files=args.input,
                out_file=args.output,
                exclude=args.exclude,
                include=args.include,
                alpha=args.alpha
            )