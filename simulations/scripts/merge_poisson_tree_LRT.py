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

    filter_muts = df['filtered'].unique()
    use_true_muts = df['true_muts'].unique()

    for filter_type in filter_muts:
        for muts_type in use_true_muts:
            total = df.shape[0]
            avg_row = np.full(df.shape[1], -1, dtype=float)

            rel_df = df[(df['filtered'] == filter_type) \
                & (df['true_muts'] == muts_type)]
            avg_row[3:-1] = rel_df.iloc[:,3:-1].mean(axis=0).values

            df.loc[total] = avg_row

            model_total = rel_df.dropna().shape[0]
            if clock:
                try:
                    avg_hyp = [f'{rel_df.iloc[:,-1].value_counts()["H0"]}/{model_total}']
                except KeyError:
                    avg_hyp = [f'0/{model_total}']
            else:
                try:
                    avg_hyp = [f'{rel_df.iloc[:,-1].value_counts()["H1"]}/{model_total}']
                except KeyError:
                    avg_hyp = [f'0/{model_total}']
            df.loc[total, 'filtered'] = int(filter_type)
            df.loc[total, 'true_muts'] = int(muts_type)
            df.loc[total, 'hypothesis_poissonTree'] = avg_hyp
    df.to_csv(out_file, sep='\t', index=False)

    avg_row = df[(df['run'] == -1) & (df['filtered'] == True) \
        & (df['true_muts'] == False)]
    print(avg_row.iloc[0])

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