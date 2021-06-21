#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_summaries(in_files, out_file):
    clock = re.search('clock(\d+.?\d*)_', out_file).group(1) == '0'
    df = pd.DataFrame()

    for i, in_file in enumerate(sorted(in_files)):
        method = os.path.basename(in_file).split('.')[0]
        if method == 'poisson':
            suffix = f'{method}'
        else:
            tree = in_file.split(os.sep)[-2].split('_')[-1]
            suffix = f'{method}.{tree}'

        new_df = pd.read_csv(in_file, sep='\t', index_col='run').dropna(axis=1)
        new_df.rename(lambda x: f'{x}_{suffix}', axis=1, inplace=True)
        df = pd.concat([df, new_df], axis=1)

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