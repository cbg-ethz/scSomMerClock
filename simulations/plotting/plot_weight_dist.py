#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd

from matplotlib.patches import Patch
from matplotlib.colors import to_rgba

from defaults import *


# colors = {101: '#FF7F00', r'$(\gamma + \beta)$ 1000': '#377DB8', 999: '#E41A1A'}

def generate_weights_plot(args):
    df = pd.read_csv(args.input, sep='\t', index_col=0)
    df = df[(df['wMax'].isin(args.wMax)) & (df['ADO'].isin(args.ADO)) \
            & (df['tree'].isin(args.method))]

    single_plot = len(args.wMax) == 1 and len(args.ADO) == 1

    wMax_vals = df['wMax'].unique()
    wMax_vals.sort()
    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig, axes = plt.subplots(nrows=wMax_vals.size, ncols=ADO_vals.size,
            figsize=(3 * ADO_vals.size, wMax_vals.size + 1))
    axes = np.reshape(axes, (wMax_vals.size, ADO_vals.size))

    for i, ADO_val in enumerate(ADO_vals):
        if single_plot:
            col_type = 'only'
        elif i == 0:
            col_type = 'first'
        elif i == ADO_vals.size - 1:
            col_type = 'last'
        elif i == np.floor(ADO_vals.size / 2):
            col_type = 'middle'
        else:
            col_type = 'intermediate'

        df_plot = df[df['ADO'] == ADO_val]
        plot_weight_dist(df_plot, wMax_vals, axes[:,i], col_type)

        if not single_plot:
            axes[0,i].annotate(f'ADO rate: {ADO_val}', xy=(0.5, 1.1),
                xytext=(0, 5), xycoords='axes fraction',
                textcoords='offset points', size='large', ha='center',
                va='baseline')

    if single_plot:
        fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.75)
    else:
        fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
            hspace=0.5, wspace=0.5)

    if args.output:
        fig.savefig(args.output, dpi=DPI)

    else:
        plt.show()
    plt.close()


def plot_weight_dist(df, wMax, axes, column_type):
    max_w = 0
    for i, j in enumerate(wMax):
        ax = axes[i]
        data = {}
        for tree, tree_data in df[df['wMax'] == j].groupby('tree'):
            w_str = tree_data['weights'].values[0].replace(';', ',')
            data[tree] = np.array(w_str.split(','), dtype=float)
            max_w = max(max_w, data[tree].max())


        sns.histplot(data,
            element='poly', stat='density', kde=False,
            common_norm=False, fill=True,
            bins=80, multiple='dodge',
            line_kws={'alpha': 0.75},
            alpha=0.4, palette=colors,
            log_scale=(0, 10),
            legend=False, ax=ax
        )

        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)
        # if i < wMax.size - 1:
        #     ax.set_xticklabels([])
        #     ax.set_xlabel(None)

        if column_type == 'first':
            if i != np.floor(wMax.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel('log Density')
        else:
            ax.set_ylabel('')

        if column_type == 'middle' and i == wMax.size - 1:
            ax.set_xlabel('Poisson Tree weights')
        else:
            ax.set_xlabel('')

        if column_type == 'last':
            ax2 = ax.twinx()
            ax2.set_ylabel('\n' + r'$w_{max}=$ '+ f'\n{j:.0f}', fontsize=12)
            ax2.set_yticks([])

        if column_type == 'only':
            ax.set_ylabel('Probability')
            ax.set_xlabel('Poisson Tree weights')

    for ax in axes:
        ax.set_xlim((0, max_w * 1.05))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input PTT weight summary file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', nargs='+', type=float,
        default=[1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
        help='wMax values to plot. Default = all.')
    parser.add_argument('-a', '--ADO', nargs='+', type=float,
        default=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
        help='ADO values to plot. Default = all.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy', 'scite'],
        default=['cellcoal', 'scite', 'cellphy'],
        help='Method to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_weights_plot(args)
