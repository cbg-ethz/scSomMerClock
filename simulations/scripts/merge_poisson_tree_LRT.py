#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_LRT(in_files, out_file):
    out_str = ''
    df = pd.DataFrame()

    clock = re.search('clock(\d+.?\d*)_', out_file).group(1) == '0'

    for i, in_file in enumerate(sorted(in_files)):
        new_df = pd.read_csv(in_file, sep='\t')
        df = df.append(new_df, ignore_index=True)
    
    # from plotting import plot_test_statistic
    # plot_test_statistic(df)

    total = df.dropna().shape[0]
    avg_row = [-1]

    for model_idx in range((len(df.columns) - 1) // 6):
        model_df = df[df.columns[6 * model_idx + 1: 6 * model_idx + 7]]
        model_df = model_df[~model_df[model_df.columns[-2]] \
            .apply(lambda x: x if isinstance(x, float) else None).isna()]
        avg_row += model_df.mean(axis=0).tolist()

        if clock:
            try:
                avg_row += [f'{model_df.iloc[:,-1].value_counts()["H0"]}/{total}']
            except KeyError:
                avg_row += [f'0/{total}']
        else:
            try:
                avg_row += [f'{model_df.iloc[:,-1].value_counts()["H1"]}/{total}']
            except KeyError:
                avg_row += [f'0/{total}']
    df.loc[total] = avg_row
    df.to_csv(out_file, sep='\t', index=False)

    print(df.loc[total])

    from plotting import plot_test_statistic, generate_pval_plot
    plot_test_statistic(df)
    generate_pval_plot(df)


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