#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_bulk_summaries(in_files, out_file):
    clock = '_clock0_' in out_file
    cols = ['R^2_pVal', 'area_pVal', 's_Bayes', 'clones_Bayes']
    if clock:
        cols.append('aff. cells')
    df = pd.DataFrame([], columns=cols)

    for i, in_file in enumerate(sorted(in_files)):
        run_no = f'{i+1:0>4}'
        with open(in_file, 'r') as f:
            log_lines = f.read().strip().split('\n')[1:]

            for j, line_raw in enumerate(log_lines):
                line = line_raw.strip()
                if line.startswith('mu'):
                    header = re.split('\s+', line)
                    if len(header) == 2: # only mu & exponent: no clone
                        df.loc[i, 's_Bayes'] = 0
                        df.loc[i, 'clones_Bayes'] = 0
                    elif 's' in header:
                        s_idx = header.index('s') + 1
                        s1 = float(re.split('\s+', log_lines[j + 1])[s_idx])
                        if log_lines[j + 2].strip().startswith('2'):
                            freq_idx = header.index('subclonefrequency') + 1
                            freq1 = float(
                                re.split('\s+', log_lines[j + 1])[freq_idx])
                            freq2 = float(
                                re.split('\s+', log_lines[j + 2])[freq_idx])
                            s2 = float(re.split('\s+', log_lines[j + 2])[s_idx])
                            df.loc[i, 's_Bayes'] = (s1 * freq1 + s2 * freq2) \
                                / (freq1 + freq2)
                            df.loc[i, 'clones_Bayes'] = 2
                        else:
                            df.loc[i, 's_Bayes'] = s1
                            df.loc[i, 'clones_Bayes'] = 1

                elif line.startswith('s'):
                    if line.endswith('subclone'):
                        freq1 = float(re.split('\s+', log_lines[j - 2])[5])
                        freq2 = float(re.split('\s+', log_lines[j - 1])[5])
                        s1 = float(re.split('\s+', log_lines[j + 1])[1])
                        s2 = float(re.split('\s+', log_lines[j + 2])[1])
                        # weighted mean
                        df.loc[i, 's_Bayes'] = (s1 * freq1 + s2 * freq2) \
                            / (freq1 + freq2)
                        df.loc[i, 'clones_Bayes'] = 2
                    else:
                        df.loc[i, 's_Bayes'] = float(
                            log_lines[j + 1].strip().split(' ')[1])
                        df.loc[i, 'clones_Bayes'] = 1
                # Best statistic, accrording to Williams et al. 2018, p. 11
                elif line == 'Area:':
                    df.loc[i, 'area_pVal'] = float(
                        log_lines[j + 1].strip().split(' ')[-1])
                elif line == 'R^2:':
                    df.loc[i, 'R^2_pVal'] = float(
                        log_lines[j + 1].strip().split(' ')[-1])
                    # Last relevant line
                    break

    df.loc[-1, 's_Bayes'] = df['s_Bayes'].mean().round(2)
    df.loc[-1, 'clones_Bayes'] = df['clones_Bayes'].mean().round(2)
    if clock:
        r2_true = (df['R^2_pVal'] >= 0.05).sum()
        area_true = (df['area_pVal'] >= 0.05).sum()
        df.loc[-1, 'aff. cells'] = -1
    else:
        r2_true = (df['R^2_pVal'] < 0.05).sum()
        area_true = (df['area_pVal'] < 0.05).sum()
        cellcoal_log = os.sep.join(
            in_file.split(os.path.sep)[:-2] + ['cellcoal.out'])
        if os.path.exists(cellcoal_log):
            with open(cellcoal_log, 'r') as f:
                for line in f:
                    if line.startswith('Data set'):
                        run_no = int(line.strip().split(' ')[2]) - 1
                    elif line.startswith('  node '):
                        cells_raw = line.strip().split('(')[-1].rstrip(')').strip()
                        cell_no = len(cells_raw.split(' '))
                        df.loc[run_no, 'aff. cells'] = cell_no
            df.loc[-1, 'aff. cells'] = df['aff. cells'].mean()
        else:
            df.loc[-1, 'aff. cells'] = -1

    df.loc[-1, 'R^2_pVal'] = f'{r2_true}/{(df["R^2_pVal"].notna()).sum()}'
    df.loc[-1, 'area_pVal'] = f'{area_true}/{(df["area_pVal"].notna()).sum()}'

    df.index.name = 'run'
    df.to_csv(out_file, sep='\t', index=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_bulk_summaries(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        merge_bulk_summaries(args.input, args.output)