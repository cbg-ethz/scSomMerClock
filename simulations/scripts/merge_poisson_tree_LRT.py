#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd


def merge_LRT_weight(in_files, out_file):
    df = pd.DataFrame()

    clock = '_clock0_' in out_file

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


def merge_LRT_tree(in_files, out_file):
    file_map = [(int(i.split('.')[-2].replace('wMax', '')), i) for i in in_files]

    for i, (_, in_file) in enumerate(sorted(file_map)):
        new_df = pd.read_csv(in_file, sep='\t', index_col=0)
        if i == 0:
            df = new_df
        else:
            df = df.join(new_df)
    df.to_csv(out_file, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-t', '--tree', action='store_true',
        help='If set, merge PoissonTree out_files over trees.'
            ' If not, merge PoissonTree out_files over weights.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if 'wMax' in snakemake.output[0]:
            merge_LRT_weight(snakemake.input, snakemake.output[0])
        else:
            merge_LRT_tree(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        if args.tree:
            merge_LRT_tree(args.input, args.output)
        else:
            merge_LRT_weight(args.input, args.output)