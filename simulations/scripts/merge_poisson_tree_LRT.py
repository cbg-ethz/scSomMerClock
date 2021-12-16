#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd


def merge_LRT(in_files, out_file):
    df = pd.DataFrame()

    clock = re.search('clock(\d+.?\d*)_', out_file).group(1) == '0'

    for i, in_file in enumerate(sorted(in_files)):
        new_df = pd.read_csv(in_file, sep='\t')
        df = df.append(new_df, ignore_index=True)

    total = df.shape[0]
    avg_row = np.full(df.shape[1], -1, dtype=float)

    avg_row[1:-2] = df.iloc[:,1:-2].mean(axis=0).values
    avg_row[-1] = np.nan

    model_total = df.dropna().shape[0]
    df.loc[total] = avg_row

    if clock:
        try:
            avg_hyp = [f'{df.iloc[:,-2].value_counts()["H0"]}/{model_total}']
        except KeyError:
            avg_hyp = [f'0/{model_total}']
    else:
        try:
            avg_hyp = [f'{df.iloc[:,-2].value_counts()["H1"]}/{model_total}']
        except KeyError:
            avg_hyp = [f'0/{model_total}']
    df.iloc[total, -2] = avg_hyp

    df.round(4).to_csv(out_file, sep='\t', index=False)

    print(df.iloc[total])

    # try:
    #     from plotting import plot_test_statistic, generate_pval_plot
    #     plot_test_statistic(df)
    #     generate_pval_plot(df)
    # except:
    #     pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_LRT(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        merge_LRT(args.input, args.output)