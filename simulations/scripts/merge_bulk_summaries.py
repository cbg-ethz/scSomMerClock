#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_bulk_summaries(in_files, out_file):
    clock = '_clock0_' in out_file
    cols = ['R^2_pVal', 'area_pVal', 's_Bayes']
    df = pd.DataFrame([], columns=cols)

    for i, in_file in enumerate(sorted(in_files)):
        run_no = f'{i+1:0>4}'
        with open(in_file, 'r') as f:
            log_lines = f.read().strip().split('\n')
            for j, line_raw in enumerate(log_lines):
                line = line_raw.strip()
                if line == 's':
                    df.loc[i, 's_Bayes'] = float(
                        log_lines[j + 1].strip().split(' ')[1])
                # Best statistic, accrording to Williams et al. 2018, p. 11
                elif line == 'Area:':
                    df.loc[i, 'area_pVal'] = float(
                        log_lines[j + 1].strip().split(' ')[-1])
                elif line == 'R^2:':
                    df.loc[i, 'R^2_pVal'] = float(
                        log_lines[j + 1].strip().split(' ')[-1])

    df.index.name = 'run'

    for col in cols:
        if col == 's_Bayes':
            df.loc[-1, col] = df[col].mean().round(2)
        else:
            if clock:
                df.loc[-1, col] = f'{(df[col] >= 0.05).sum()}/{(df[col].notna()).sum()}'
            else:
                df.loc[-1, col] = f'{(df[col] < 0.05).sum()}/{(df[col].notna()).sum()}'
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