#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_LRT(in_files, out_file):
    out_str = ''
    df = pd.DataFrame()

    clock = re.search('clock(\d+.?\d*)_', out_file).group(1) == '0'

    for i, in_file in enumerate(in_files):
        new_df = pd.read_csv(in_file, sep='\t')
        df = df.append(new_df, ignore_index=True)
    
    means = df.mean(axis=0)
    total = df.shape[0]
    if clock:
        try:
            H_poisson = f'{df["hypothesis_poisson"].value_counts()["H0"]}/{total}'
        except KeyError:
            H_poisson = f'0/{total}'
        try:
            H_binom = f'{df["hypothesis_nbinom"].value_counts()["H0"]}/{total}'
        except KeyError:
            H_binom = f'0/{total}'
    else:
        try:
            H_poisson = f'{df["hypothesis_poisson"].value_counts()["H1"]}/{total}'
        except KeyError:
            H_poisson = f'0/{total}'
        try:
            H_binom = f'{df["hypothesis_nbinom"].value_counts()["H1"]}/{total}'
        except KeyError:
            H_binom = f'0/{total}'
    
    avg_row = [-1] + means.iloc[1:6].tolist() + [H_poisson] \
        + means.iloc[6:].tolist() + [H_binom]
    df.loc[total] = avg_row

    df.to_csv(out_file, sep='\t', index=False)

    print(df.loc[total])


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