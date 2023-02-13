#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *
AFF_MIN = 0.1
AFF_MAX = 0.9


def generate_pval_plot_ss(args):
    cols =  ['amplifier', 'ss_size', 'ss_rep', 'ss_aff', 'method', 'P-value']
    df = pd.DataFrame(columns=cols)

    for res_file in sorted(os.listdir(args.input)):
        if not 'ss_summary' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO != args.ADO:
            continue

        ampl = float(re.search('res_clock(\d[\.\d]*)', res_file).group(1))
        if ampl == 0:
            ampl = 1

        if ampl not in args.amplifier:
            continue

        print(f'Including summary file: {res_file}')

        df_summary = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_summary.drop([-1], inplace=True)

        rel_cols = ['subsample_size', 'subsample_rep']
        if ampl > 1:
            df_summary = df_summary[
                (df_summary['aff. cells'] >= args.total_cells * args.aff_min) \
                & (df_summary['aff. cells'] <= args.total_cells * args.aff_max)]

            rel_cols.append('aff. cells sampled')

        for col in df_summary.columns:
            if col.startswith('p-value'):
                wMax = int(col.split('_')[-2].replace('wMax', ''))
                if wMax != args.wMax:
                    continue

                tree = col.split('_')[-1]
                if tree not in args.method:
                    continue

                new_df = df_summary[rel_cols + [col]]
                new_df.insert(0, 'method', tree)
                new_df.insert(0, 'amplifier', ampl)
                new_df = new_df.rename({col: 'P-value', 'subsample_size': 'ss_size',
                    'subsample_rep': 'ss_rep', 'aff. cells sampled': 'ss_aff'},
                    axis=1)

                df = df.append(new_df, ignore_index=True)

    df = df[df['ss_size'].isin(args.subsamples)]

    subsamples = df['ss_size'].unique()
    subsamples.sort()
    ampl_vals = df['amplifier'].unique()
    ampl_vals.sort()

    col_no = subsamples.size
    row_no = ampl_vals.size
    fig, axes = get_subplots(row_no, col_no)

    for i, ampl in enumerate(ampl_vals):
        if ampl == 1:
            row_title = 'Clock'
        else:
            row_title = f'Rate change:\n{ampl:.0f}x'

        df_plot = df[(df['amplifier'] == ampl)]
        plot_pVal_dist(df_plot, subsamples, axes[i], (i, row_no), row_title)

    if col_no > 1:
        for i, ss in enumerate(subsamples):
            add_col_header(axes[0, i], f'Cells:\n{ss}/{args.total_cells}')

    plot_fig(fig, args.output)


def plot_pVal_dist(df, subsamples, axes, row_no, row_title=''):
    hue_order = ['cellcoal', 'cellphy']
    for col_no, ss in enumerate(subsamples):
        ax = axes[col_no]
        data = df[df['ss_size'] == ss]
        if row_title.lower() == 'clock':
            n_annot = r'$n$' + f' = {data["method"].value_counts().min()}'
        else:
            no_ampl = np.sum(data['ss_aff'] == 0) / data["method"].unique().size
            data = data[data['ss_aff'] > 0]
            if len(str(no_ampl)) == 1:
                n_annot = r'$n$' + f'  = {data["method"].value_counts().min(): >3}\n' \
                    + r'$n_0$' + f' = {no_ampl: >5.0f}'
            else:
                n_annot = r'$n$' + f'  = {data["method"].value_counts().min(): >3}\n' \
                    + r'$n_0$' + f' = {no_ampl: >4.0f}'

        dp = sns.histplot(data, x='P-value', hue='method', ax=ax, **HIST_DEFAULT)

        ax.set_xlim((0, 1))
        ax.set_xticks([0.0, 0.5, 1])
        ax.set_ylim((0, 1))
        ax.set_yticks([0.0, 0.5, 1])

        # Add rugplots
        k = 0
        for method in hue_order:
            rug_data = data[data['method'] == method]['P-value'].values
            if rug_data.size == 0:
                continue
            add_rugs(rug_data, offset=k, ax=ax, color=COLORS[method])
            k += 1

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.annotate(n_annot, xy=(0.9, 0.825),
            xytext=(0, 5), xycoords='axes fraction', textcoords='offset points',
            ha='right', va='top', bbox=bbox_props)


        # Only x axis on last row
        if row_no[0] < row_no[1] - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        # y label on first col, mid row
        if col_no == 0:
            if row_no[0] != np.floor(row_no[1] / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel('Probability')
        else:
            ax.set_ylabel(None)
            ax.set_yticklabels([])

        # x label on mid col, last row,
        if col_no == np.floor(subsamples.size / 2) and row_no[0] == row_no[1] - 1:
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')

        if col_no == subsamples.size - 1:
            add_row_header(ax, f'\n{row_title}')


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with subsampling')
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', type=float, default=400,
        help='wMax value to plot. Default = 400.')
    parser.add_argument('-a', '--amplifier', nargs='+', default=[1, 2, 5],
        type=int, help='Amplifier values to plot. Default = all.')
    parser.add_argument('--ADO', default=0.2, type=float,
        help='Simulated ADO value to plot. Default = 0.2.')
    parser.add_argument('-ss', '--subsamples', nargs='+', type=int,
        choices=[10, 30, 50, 70, 90], default=[10, 30, 50, 70, 90],
        help='# Cells subsampled. Default = [10, 30, 50, 70, 90].')
    parser.add_argument('-t', '--total_cells', type=int, default=100,
        help='Total number of simualted cells. Default = 100.')
    parser.add_argument('--aff_min', type=float, default=0.1,
        help='Min. number of affected cells for a run to be considered. ' \
            'Default = 0.1.')
    parser.add_argument('--aff_max', type=float, default=0.9,
        help='Max. number of affected cells for a run to be considered. ' \
            'Default = 0.9.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy'], default=['cellcoal','cellphy'],
        help='Method to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_ss(args)