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
        subs_re = re.search('\.ss(\d+)\.', in_file)
        if subs_re:
            new_df['subsample_size'] = int(subs_re.group(1))
            new_df['subsample_rep'] = int(
                re.search('\.ss\d+\.(\d+)\.wMax', in_file).group(1))
        df = df.append(new_df, ignore_index=True)

    num_cols = df.columns[df.dtypes != object]
    mean_col = df.loc[:,num_cols].mean(axis=0)

    model_total = df.dropna().shape[0]
    hyp_col = [i for i in df.columns if i.startswith('hypothesis')][0]
    if clock:
        try:
            avg_hyp = f'{df[hyp_col].value_counts()["H0"]}/{model_total}'
        except KeyError:
            avg_hyp = f'0/{model_total}'
    else:
        try:
            avg_hyp = f'{df[hyp_col].value_counts()["H1"]}/{model_total}'
        except KeyError:
            avg_hyp = f'0/{model_total}'

    df.loc[-1] = np.full(df.shape[1], -1, dtype=float)
    df.loc[-1, num_cols] = mean_col
    df.loc[-1, 'run'] = -1
    if subs_re:
        df.loc[-1, 'subsample_size'] = ','.join(
            [f'{i:.0f}' for i in df['subsample_size'].unique()])
        df.loc[-1, 'subsample_rep'] = df['subsample_rep'].max()

    df.loc[-1, hyp_col] = avg_hyp
    df.round(4).to_csv(out_file, sep='\t', index=False)
    print(df.loc[-1])


def merge_LRT_tree(in_files, out_file):
    df = pd.DataFrame()
    clock = '_clock0_' in out_file

    for in_file in sorted(in_files):
        new_df = pd.read_csv(in_file, sep='\t')
        # Backward compatibility
        if 'SNVs' in new_df.columns:
            new_df.drop(['SNVs', 'TP', 'FP', 'TN', 'FN', 'MS', 'MS_T'],
                axis=1, errors='ignore', inplace=True)

        subs_re = re.search('\.ss(\d+)\.', in_file)
        if subs_re:
            new_df['subsample_size'] = int(subs_re.group(1))
            new_df['subsample_rep'] = int(
                re.search('\.ss\d+\.(\d+)\.LRT', in_file).group(1))

        df = df.append(new_df, ignore_index=True)

    if 'subsample_size' in df.columns:
        df.reset_index(inplace=True, drop=True)
        df.set_index(['run', 'subsample_size', 'subsample_rep'], inplace=True)
        idx_row = (-1, -1, -1)
    else:
        df.set_index('run', inplace=True, drop=True)
        idx_row = -1

    df.sort_index(inplace=True)
    # Add average row
    num_cols = [i for i, (j, k) in enumerate(df.dtypes.items()) if k != object]
    mean_col = df.iloc[:,num_cols].mean(axis=0)

    model_total = df.dropna().shape[0]
    if isinstance(df.index, pd.MultiIndex):
        avg_id = (-1, -1, -1)
    else:
        avg_id = -1
    df.loc[avg_id] = np.full(df.shape[1], -1, dtype=float)
    df.loc[avg_id] = mean_col

    for i, hyp_col in enumerate(df.columns):
        if not hyp_col.startswith('hypothesis'):
            continue

        if clock:
            try:
                avg_hyp = f'{df[hyp_col].value_counts()["H0"]}/{model_total}'
            except KeyError:
                avg_hyp = f'0/{model_total}'
        else:
            try:
                avg_hyp = f'{df[hyp_col].value_counts()["H1"]}/{model_total}'
            except KeyError:
                avg_hyp = f'0/{model_total}'

        df.loc[avg_id, i] = avg_hyp

    df.round(4).to_csv(out_file, sep='\t')
    print(df.loc[avg_id])


def merge_LRT_weight_column(in_dir, out_file=''):
    vals = []
    ADO = float(re.search('WGA(0[\.\d]*)', in_dir).group(1))

    for tree in ['cellcoal', 'cellphy', 'scite']:
        in_file = os.path.join(in_dir, f'poissonTree_{tree}',
            'poissonTree.summary.tsv')
        df_in = pd.read_csv(in_file, sep='\t', index_col=0)
        df_in.drop([-1], inplace=True)
        for col in df_in.columns:
            if not col.startswith('weights'):
                continue
            wMax = float(col.split('_')[-1].replace('wMax', ''))
            if wMax == 1:
                continue
            vals.append([ADO, tree, wMax, ';'.join(df_in[col].values)])

    df = pd.DataFrame(vals, columns=['ADO', 'tree', 'wMax', 'weights'])

    if not out_file:
        out_file = os.path.join(in_dir, 'PTT_weights.tsv')
    df.to_csv(out_file, sep='\t')


def get_LRT_weight_summary(in_dir, out_file=''):
    df = pd.DataFrame(columns=['ADO', 'tree', 'wMax', 'weights'])
    for dir in os.listdir(in_dir):
        if not dir.startswith('res_clock0') or 'bulk' in dir:
            continue
        in_file = os.path.join(in_dir, dir, 'minDP5-minGQ1', 'PTT_weights.tsv')
        if not os.path.exists(in_file):
            print(f'!Warning: Missing file {in_file}')
            continue
        new_df = pd.read_csv(in_file, sep='\t', index_col=0)
        df = df.append(new_df, ignore_index=True)

    if not out_file:
        out_file = os.path.join(in_dir, 'PTT_weights_all.tsv')
    df.to_csv(out_file, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-t', '--tree', action='store_true',
        help='If set, merge PoissonTree out_files over trees. '
            'If not, merge PoissonTree out_files over weights.')
    parser.add_argument('-w', '--weight', action='store_true',
        help='Generate weight file for single ADO rate.')
    parser.add_argument('-ws', '--weight_summary', action='store_true',
        help='Generate weight summary file over all ADO rates.')
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

        if args.weight:
            merge_LRT_weight_column(args.input[0], args.output)
            exit()

        if args.weight_summary:
            get_LRT_weight_summary(args.input[0], args.output)
            exit()

        if args.tree:
            merge_LRT_tree(args.input, args.output)
        else:
            merge_LRT_weight(args.input, args.output)