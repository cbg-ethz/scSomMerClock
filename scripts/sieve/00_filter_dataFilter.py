#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
from itertools import cycle
from matplotlib import pyplot as plt


def main(in_file, exome_file, out_file):
    df = pd.read_csv(exome_file, sep='\t', names=['chr', 'start', 'stop'],
        dtype={'chr': str})
    df['chr'].replace({'X': '23', 'Y': '24'}, inplace=True)
    df = df[df['chr'].apply(lambda x: x.isnumeric())]
    df['chr'] = df['chr'].astype(int)
    df.sort_values(['chr', 'start', 'stop'], axis=0, ignore_index=True,
        inplace=True)
    import pdb; pdb.set_trace()

    df_sum.to_csv(out_file, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='SIEVE output file')
    parser.add_argument('exome', type=str, help='Exome region file (bed format)')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output File. Default = <INPUT>.filtered.tsv'
    )
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_summaries(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = f'{os.path.splitext(args.input)[0]}.filtered.tsv'
        main(args.input, args.exome, args.output)