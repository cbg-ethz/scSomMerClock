#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *


def generate_power_plot(args):
    cols =  ['ADO', 'amplifier', 'ss_size', 'method', 'power']
    df = pd.DataFrame(columns=cols)

    for res_file in sorted(os.listdir(args.input)):
        if not os.path.isfile(os.path.join(args.input, res_file)):
            continue
        if not res_file.startswith('res_clock') or not 'bulk' in res_file:
            continue

        if '_ss_' in res_file:
            ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
            if ADO not in args.ADO:
                continue
        else:
            ADO = 0

        ampl = float(re.search('res_clock(\d[\.\d]*)', res_file).group(1))
        if ampl not in args.amplifier:
            continue

        print(f'Including summary file: {res_file}')

        try:
            df_new = pd.read_csv(
                os.path.join(args.input, res_file), sep='\t', index_col='run')
        except:
            import pdb; pdb.set_trace()
        df_new.drop([-1], inplace=True)

        if ampl > 1:
            if '_ss_' in res_file:
                min_frac = args.clone_size[0]
                max_frac = args.clone_size[1]
            else:
                min_frac = 0.2
                max_frac = 0.7
            df_new = df_new[(df_new['aff. cells'] >= args.total_cells * min_frac) \
                & (df_new['aff. cells'] <= args.total_cells * max_frac)]

        if '_ss_' in res_file:
            for ss_size, ss_df in df_new.groupby('subsample_size'):
                if ss_size not in args.subsamples:
                    continue

                for col in ss_df.columns:
                    if not col.startswith('hypothesis'):
                        continue

                    tree = re.split('[\._]', col)[-1]
                    if tree not in args.method:
                        continue

                    wMax = int(re.search('wMax(\d+)', col).group(1))
                    if wMax != args.wMax:
                        continue

                    power = (ss_df[col] == 'H1').mean()
                    df.loc[df.shape[0]] = [ADO, ampl, ss_size, tree, power]
        else:
            power = (df_new['area_pVal'].astype(float) < 0.05).mean()
            df.loc[df.shape[0]] = [ADO, ampl, args.total_cells,
                    'neutrality', power]

    row_no = 1
    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()
    col_no = ADO_vals.size

    fig, axes = get_subplots(row_no, col_no)

    for i, ADO_val in enumerate(ADO_vals):
        df_plot = df[df['ADO'] == ADO_val]
        plot_power(df_plot, axes[0, i], args.total_cells, (i, col_no))
        if col_no > 1:
            if ADO_val > 0:
                col_title = f'FN rate: {ADO_val / 2}'
            else:
                col_title = 'No Errors'
            add_col_header(axes[0, i], col_title)

    plot_fig(fig, args.output)


def plot_power(df, ax, n_total, col):
    labels = []
    handles = []

    for method, method_df in df.groupby('method'):
        for ss, ss_df in method_df.groupby('ss_size'):
            x = ss_df['amplifier'].values
            y = ss_df['power'].values
            labels.append(f'{method} ({ss} / {n_total} cells)')
            print(method, '\n', ss_df, '\n\n')
            sns.lineplot(x=x, y=y, ax=ax, color=COLORS[method],
                linestyle=LINE_STYLE[ss], linewidth=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # First col
    if col[0] == 0:
        ax.set_ylabel(r'$P(H_0) \leq 0.05$ [%]')
        ax.set_yticklabels([0, 50, 100])
    else:
        ax.set_ylabel(None)
        ax.set_yticklabels([])

    ax.set_yticks([0, .50, 1])
    ax.set_ylim((0, 1.05))

    # Middle col
    if col[0] == np.floor(col[1] / 2):
        ax.set_xlabel('Rate change')
    else:
        ax.set_xlabel('')

    ax.set_xticks(df['amplifier'].unique())
    ax.set_xticklabels([f'{i: >2.0f}x' if i > 1 else 'Clock' \
        for i in df['amplifier'].unique()], rotation=90)


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with branch multiplier')
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', type=float, default=500,
        help='wMax value to plot. Default = 500.')
    parser.add_argument('-a', '--ADO', nargs='+', type=float,
        default=DEF_ADO,
        help=f'ADO values to plot. Default = {DEF_ADO}.')
    parser.add_argument('-amp', '--amplifier', nargs='+', default=[0, 2, 5, 10],
        type=float, help='Amplifier value to plot. Default = [1, 2, 5, 10]')
    parser.add_argument('-c', '--clone_size', default = [0.1, 0.9],
        type=float, help='Amplified clone size subsets. Default = [0.1, 0.9].')
    parser.add_argument('-ss', '--subsamples', nargs='+', type=int,
        choices=[10, 30, 50, 70, 90], default=[10, 50, 90],
        help='# Cells subsampled. Default = [10, 50, 90].')
    parser.add_argument('-t', '--total_cells', type=int, default=100,
        help='Total number of simulated cells. Default = 100.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy', 'scite', 'neutralitytest'],
        default=['cellcoal', 'cellphy', 'mobster', 'neutralitytest'],
        help='Method to plot.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_power_plot(args)