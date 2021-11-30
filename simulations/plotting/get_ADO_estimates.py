#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import re


def collect_ADO_estimates(in_dir, out_file):
    cols = ['CellPhy_FP', 'CellPhy_FN', 'Scite_FP', 'Scite_FN']
    df = pd.DataFrame([], index=[], columns=cols)

    scite_dir = os.path.join(in_dir, 'scite_dir')
    if os.path.exists(scite_dir):
        for scite_file in os.listdir(scite_dir):
            if not scite_file.endswith('.log'):
                continue

            with open(os.path.join(scite_dir, scite_file), 'r') as f:
                log = f.read()
            FN = float(re.search('best value for beta:\\\\t(\d.\d+(e-\d+)?)', log) \
                .group(1))
            FP = float(re.search('best value for alpha:\\\\t(\d.\d+(e-\d+)?)', log) \
                .group(1))
            run_no = int(scite_file.split('.')[1])
            df.loc[run_no, 'Scite_FP'] = FP
            df.loc[run_no, 'Scite_FN'] = FN

    cellphy_dir = os.path.join(in_dir, 'cellphy_dir')
    if os.path.exists(cellphy_dir):
        for cellphy_file in os.listdir(cellphy_dir):
            if not cellphy_file.endswith('.log'):
                continue

            with open(os.path.join(cellphy_dir, cellphy_file), 'r') as f:
                    log = f.read().strip()
            FP = float(re.search('SEQ_ERROR: (0.\d+(e-\d+)?)', log).group(1))
            FN = float(re.search('ADO_RATE: (0.\d+(e-\d+)?)', log).group(1))

            run_no = int(cellphy_file.split('.')[1])
            df.loc[run_no, 'CellPhy_FP'] = FP
            df.loc[run_no, 'CellPhy_FN'] = FN

    df.index.name = 'run'
    df.dropna(axis=1).to_csv(out_file, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if not args.output:
        args.output = os.path.join(args.input, 'ADO_estimates.tsv')
    collect_ADO_estimates(args.input, args.output)