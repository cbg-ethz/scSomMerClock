#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *


def generate_pval_plot_ss(args):
    cols =  ['ADO', 'ss_size', 'ss_rep', 'method', 'P-value']
    df = pd.DataFrame(columns=cols)

    for res_file in sorted(os.listdir(args.input)):
        if not 'ss_summary' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO not in args.ADO:
            continue

        ampl = float(re.search('res_clock(\d[\.\d]*)', res_file).group(1))
        if ampl != args.amplifier:
            continue

        print(f'Including summary file: {res_file}')

        df_summary = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_summary.drop([-1], inplace=True)

        if ampl > 0:
            df_summary = df_summary[
                (df_summary['aff. cells'] >= args.total_cells * 0.1) \
                & (df_summary['aff. cells'] <= args.total_cells * 0.9)]

        rel_cols = ['subsample_size', 'subsample_rep']
        for col in df_summary.columns:
            if col.startswith('p-value'):
                wMax = int(col.split('_')[-2].replace('wMax', ''))
                if wMax != args.wMax:
                    continue

                tree = col.split('_')[-1]
                if tree not in args.method:
                    continue

                new_df = df_summary[rel_cols + [col]]
                new_df.loc[:,'ADO'] = ADO
                new_df.loc[:,'method'] = tree
                new_df.rename({col: 'P-value', 'subsample_size': 'ss_size',
                    'subsample_rep': 'ss_rep'}, axis=1, inplace=True)

                df = df.append(new_df, ignore_index=True)

    df = df[df['ss_size'].isin(args.subsamples)]

    subsamples = df['ss_size'].unique()
    subsamples.sort()
    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()

    single_plot = subsamples.size == 1 and ADO_vals.size == 1

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(2, 2))
    else:
        fig, axes = plt.subplots(nrows=subsamples.size, ncols=ADO_vals.size,
            figsize=(2 * ADO_vals.size, subsamples.size + 1))
    axes = np.reshape(axes, (subsamples.size, ADO_vals.size))

    for i, ss in enumerate(subsamples):
        df_plot = df[(df['ss_size'] == ss)]
        plot_pVal_dist(df_plot, ADO_vals, axes[i], (i, subsamples.size),
            f'{ss}/{args.total_cells}\ncells')

    if not single_plot:
        for i, ADO_val in enumerate(ADO_vals):
            if ADO_val == 0:
                col_title = 'No Errors'
            else:
                col_title = f'FN rate: {ADO_val / 2}'

            axes[0, i].annotate(col_title, xy=(0.5, 1.1), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_pVal_dist(df, ADOs, axes, row_no, row_title=''):
    hue_order = ['cellcoal', 'cellphy']
    for col_no, j in enumerate(ADOs):
        ax = axes[col_no]
        data = df[df['ADO'] == j]
        dp = sns.histplot(data, x='P-value', hue='method', element='bars',
            stat='probability', multiple='layer', fill=True, common_norm=False,
            hue_order=hue_order, binwidth=0.05,
            binrange=(0, 1), palette=colors, shrink=1, legend=False, ax=ax,
        )

        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))

        # Add rugplots
        k = 0
        for method in hue_order:
            rug_data = data[data['method'] == method]['P-value'].values
            if rug_data.size == 0:
                continue
            add_rugs(rug_data, offset=k, ax=ax, color=colors[method])
            k += 1

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

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
        if col_no == np.floor(ADOs.size / 2)  and row_no[0] == row_no[1] - 1:
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')

        if col_no == ADOs.size - 1:
            ax2 = ax.twinx()
            ax2.set_ylabel(f'\n{row_title}', fontsize=12)
            ax2.set_yticks([])
            for tick in  ax.yaxis.majorTicks:
                tick.tick1line.set_markersize(0)


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with subsampling')
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', type=float, default=400,
        help='wMax values to plot. Default = all.')
    parser.add_argument('-a', '--amplifier', default=2, type=int,
        help='Amplifier value to plot. Default = 2.')
    parser.add_argument('-do', '--ADO', default=[0, 0.2, 0.4], type=float,
        help='Simulated ADO value to plot. Default = [0, 0.2, 0.4].')
    parser.add_argument('-ss', '--subsamples', nargs='+', type=int,
        choices=[10, 30, 50, 70, 90], default=[10, 30, 50, 70, 90],
        help='Min. #cells affected. Default = all.')
    parser.add_argument('-t', '--total_cells', type=int, default=100,
        help='Total number of simualted cells. Default = 100.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy'], default=['cellcoal','cellphy'],
        help='Method to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_ss(args)