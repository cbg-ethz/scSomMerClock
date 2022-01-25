#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_summaries(in_files, out_file):
    clock = '_clock0_' in out_file
    df = pd.DataFrame()

    for i, in_file in enumerate(sorted(in_files)):
        method = os.path.basename(in_file).split('.')[0]
        if method == 'poissonDisp':
            suffix = f'{method}'
        else:
            tree_dir = in_file.split(os.sep)[-2]
            if 'cellcoal' in tree_dir:
                tree = 'cellcoal'
            elif 'cellphy' in tree_dir:
                tree = 'cellphy'

            elif 'scite' in tree_dir:
                tree = 'scite'
            else:
                raise RuntimeError(f'Cannot determine tree from dir: {tree_dir}')

            if tree_dir.count('_') == 2:
                weight = tree_dir.split('_')[1]
                suffix = f'{method}_{weight}.{tree}'
            else:
                suffix = f'{method}.{tree}'

        new_df = pd.read_csv(in_file, sep='\t', index_col='run').dropna(axis=1)
        if method == 'poissonTree':
            to_drop = [i for i in new_df.columns \
                if 'nlopt' in i or 'multinomial' in i]
            new_df.drop(to_drop, inplace=True, axis=1)
            new_df.rename(lambda x: x.split('_')[0], axis=1, inplace=True)
        new_df.rename(index={'Avg.': '-1'}, inplace=True)
        new_df.index = new_df.index.astype(int)
        new_df.rename(lambda x: f'{x.split(":")[0]}_{suffix}', axis=1, inplace=True)

        df = pd.concat([df, new_df], axis=1)

    if not clock:
        cellcoal_log = os.path.join(in_file.split(os.path.sep)[0], 'cellcoal.out')
        with open(cellcoal_log, 'r') as f:
            for line in f:
                if line.startswith('Data set'):
                    run_no = int(line.strip().split(' ')[2])
                elif line.startswith('  node '):
                    cells_raw = line.strip().split('(')[-1].rstrip(')').strip()
                    cell_no = len(cells_raw.split(' '))
                    df.loc[run_no, 'aff. cells'] = cell_no
        df.loc[-1, 'aff. cells'] = df['aff. cells'].mean()

    snv_cols = [i for i in df.columns if i.startswith('SNVs_')]
    dof_cols = [i for i in df.columns if i.startswith('dof_')]
    col_drop = [i for i in df.columns \
        if i.startswith('H0_') or i.startswith('H1_')]

    df['dof'] = df[dof_cols[0]]
    df['#SNV'] = df[snv_cols[0]]
    df.drop(snv_cols + col_drop + dof_cols, inplace=True, axis=1)
    if not clock:
        df = df[['aff. cells', 'dof', '#SNV'] + list(df.columns[:-3])]
    else:
        df = df[['dof', '#SNV'] + list(df.columns[:-2])]

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
        merge_summaries(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        merge_summaries(args.input, args.output)