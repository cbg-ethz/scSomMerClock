#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd

from defaults import *


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

    if not out_file:
        out_file = os.path.join(in_dir, 'ADO_estimates.tsv')

    df.dropna(axis=1).to_csv(out_file, sep='\t')


def plot_ADO_estimate(in_dir, out_file='', bulk=False):
    df = pd.DataFrame(
        columns=['True ADO', 'Tree', 'Missing', 'Inferred ADO'])
    # cmap =cm.get_cmap('inferno_r') #plasma, magma, cvidis

    for res_file in sorted(os.listdir(in_dir)):
        if  not res_file.startswith('res_clock0'):
            continue

        if any([True for i in [0.8] if f'WGA{i}' in res_file]): continue

        if bulk and not 'bulk' in res_file:
            continue
        elif not bulk and 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO == 0:
            continue

        res_df = pd.read_csv(os.path.join(in_dir, res_file), sep='\t', index_col=0)
        avg_row = res_df.loc[-1]
        res_df.drop(-1, inplace=True)

        ADO_file = os.path.join(in_dir, f'ADO_overview_WGA{ADO}.tsv')
        ADO_df = pd.read_csv(ADO_file, sep='\t', index_col=0)
        ADO_df.drop('healthycell', axis=1, inplace=True)
        ADO_mean = ADO_df.mean(axis=1)
        ADO_mean.name = 'True ADO'

        OV_file = os.path.join(in_dir, f'data_overview_WGA{ADO}.tsv')
        OV_df = pd.read_csv(OV_file, sep='\t', index_col=0, encoding='utf-8')

        MS = OV_df['MS_sides'] / OV_df['sides']
        MS.name = 'Missing'

        for col, content in res_df.items():
            if not col.startswith('FN_poissonTree'):
                continue
            tree = col.split('_')[-1]
            if tree == 'cellcoal':
                continue
            new_df = pd.DataFrame(content)
            new_df.rename({col: 'Inferred ADO'}, axis=1, inplace=True)
            new_df['Tree'] = METHODS[tree]
            new_df['Missing'] = 'not added'

            new_df = new_df.merge(ADO_mean, left_index=True, right_index=True)
            df = df.append(new_df, ignore_index=True)

            new_df['Inferred ADO'] = new_df['Inferred ADO'] + MS
            new_df['Missing'] = 'added'
            df = df.append(new_df, ignore_index=True)

    df['Inferred ADO'] *= 2
    fig, ax = plt.subplots(figsize=(6, 6))

    sns.scatterplot(data=df, x='True ADO', y='Inferred ADO', style='Missing',
        hue='Tree', alpha=0.2, palette=colors, ax=ax
    )
    plt.axline((0, 0), slope=1, color='black')

    max_x = np.ceil(df['True ADO'].max() * 10) / 10

    ax.set_ylim((-0.01, max_x))
    ax.set_xlim((-0.01, max_x))

    ax.set_ylabel('Inferred ADO [%]')
    ax.set_xlabel('True ADO [%]')

    handles, labels = ax.get_legend_handles_labels()

    ax.legend(handles, labels, ncol=2, frameon=True)

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95)
    if out_file:
        fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-b', '--bulk', action='store_true',
        help='Consider only bulk file.')
    parser.add_argument('-c', '--collect', action='store_true',
        help='Generate estimate overview file from cellcoal outdir.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.collect:
        collect_ADO_estimates(args.input, args.output)
    else:
        plot_ADO_estimate(args.input, args.output, args.bulk)