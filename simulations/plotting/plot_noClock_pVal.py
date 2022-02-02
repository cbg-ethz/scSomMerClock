#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


TICK_FONTSIZE = 12
LABEL_FONTSIZE = 16
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'font.size': TICK_FONTSIZE,
        'axes.titlesize': LABEL_FONTSIZE,
        'axes.labelsize': LABEL_FONTSIZE,
        'axes.titlesize': LABEL_FONTSIZE,
        'axes.labelticksize': LABEL_FONTSIZE,
        'lines.linewidth': 1,
        'legend.fontsize': LABEL_FONTSIZE,
        'legend.title_fontsize':  LABEL_FONTSIZE,
        'xtick.major.size':  2,
        'ytick.major.size':  2,
})

colors = {
    'cellcoal': '#E41A1A', # red
    'cellphy': '#377DB8', # blue
    'scite': '#FF7F00', # orange
    '-': '#994EA3'
}


def plot_affected_cells(df, out_file):
    x = np.arange(1, df['aff. cells'].max() + 1, 1)
    y = np.zeros(x.size)
    for i, x_i in enumerate(x):
        y[i] = df[df['aff. cells'] == x_i].shape[0]
    y /= df.shape[0]

    col1 = '#377eb8'
    fig, ax = plt.subplots(figsize=(12, 9))
    sns.lineplot(x=x, y=y, ax=ax, lw=2, marker='o', markersize=8, color=col1)

    ax.set_ylim((0, y.max() * 1.05))
    y_ticks = np.arange(0, 0.01 + np.ceil(y.max() * 10) / 10, 0.05)
    ax.set_yticks(y_ticks)
    ax.set_ylabel(f'Rate of samples [%]', fontsize=LABEL_FONTSIZE)
    ax.set_xlim((0, 30))
    ax.set_xlabel(f'#Cells affected by switch', fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
    ax.spines['left'].set_color(col1)
    ax.tick_params(axis='y', colors=col1)

    col2 = '#e41a1c'
    ax2 = plt.twinx()
    y2 = 1 - np.cumsum(y)
    sns.lineplot(x=x, y=y2, ax=ax2, lw=2, marker='o', markersize=8, color=col2)
    ax2.set_ylim((0, y2.max() * 1.05))
    y2_ticks = np.linspace(0, np.ceil(y2.max() * 10) / 10, y_ticks.size)
    ax2.set_yticks(y2_ticks)
    ax2.set_yticklabels(np.round(y2_ticks, 2))
    ax2.set_ylabel(f'1 - Cum. Sum(Rate of samples)', fontsize=LABEL_FONTSIZE)
    ax2.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax2.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
    ax2.spines['right'].set_color(col2)
    ax2.tick_params(axis='y', colors=col2)

    if out_file:
        fig.savefig(os.path.splitext(out_file)[0] + '_affectedCells.png',
            dpi=300)
    else:
        plt.show()
    plt.close()


def plot_wmax_pval(in_file, out_file='', wMax_show=[]):
    df_in = pd.read_csv(in_file, sep='\t', index_col=0)
    df_in.drop([-1], inplace=True)
    rel_cols = ['aff. cells'] + [i for i in df_in.columns if 'p-value' in i]
    df_in = df_in[rel_cols]

    # plot_affected_cells(df, out_file)
    vals = []
    for i, run_data in df_in.iterrows():
        for name, value in run_data.items():
            if name == 'aff. cells':
                cell_no = value
            elif name == 'p-value_poissonDisp':
                vals.append(['-', -1, cell_no, value])
            else:
                tree = name.split('.')[-1]
                wMax = int(re.search('_wMax(\d+)', name).group(1))
                vals.append([tree, wMax, cell_no, value])
    df = pd.DataFrame(vals,
        columns=['Tree', 'wMax', 'aff. cells', 'P-value'])

    if wMax_show:
        df = df[df['wMax'].isin(wMax_show)]

    min_cells = np.array([1, 2, 3, 4, 5])
    wMax_vals = df['wMax'].unique()
    wMax_vals.sort()

    fig, axes = plt.subplots(nrows=wMax_vals.size, ncols=min_cells.size,
        figsize=(3 * min_cells.size, wMax_vals.size + 1))

    for i, min_cell in enumerate(min_cells):
        if i == 0:
            col_type = 'first'
        elif i == min_cells.size - 1:
            col_type = 'last'
        elif i == np.floor(min_cells.size / 2):
            col_type = 'middle'
        else:
            col_type = 'intermediate'

        plot_pVal_dist(df[df['aff. cells'] >= min_cell], wMax_vals, axes[:,i],
            col_type, min_cell)

        header = r'$\geq$' + f'{min_cell} cells\naffected'
        axes[0,i].annotate(header, xy=(0.5, 1), xytext=(0, 5),
            xycoords='axes fraction', textcoords='offset points',
            size='large', ha='center', va='baseline')

    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
        hspace=0.5, wspace=0.5)

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_pVal_dist(df, wMax, axes, column_type, min_cell):
    for i, j in enumerate(wMax):
        ax = axes[i]
        dp = sns.histplot(df[df['wMax'] == j], x='P-value', hue='Tree',
            element='poly', stat='probability', kde=False,
            common_norm=False, fill=False,
            binwidth=0.01, binrange=(0, 1), multiple='dodge',
            kde_kws={'cut': 0, 'clip': (0, 1)},
            line_kws={'lw': 3, 'alpha': 0.66},
            palette=colors, alpha=0.5,
            hue_order=['-', 'cellcoal', 'cellphy', 'scite'],
            legend=False, ax=ax,
        )
        for tree, col in colors.items():
            df_sub = df[(df['Tree'] == tree) & (df['wMax'] == j)]
            if df_sub.size > 0:
                ax.axhline((df_sub['P-value'] <= 0.05).mean(), ls='--',
                    color=col, lw=1)

        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < wMax.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if column_type == 'first':
            if i != np.floor(wMax.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel(f'Probability', fontsize=LABEL_FONTSIZE)
        else:
            ax.set_ylabel('')

        if column_type == 'middle' and i == wMax.size - 1:
            ax.set_xlabel('P-value', fontsize=LABEL_FONTSIZE)
        else:
            ax.set_xlabel('')

        if column_type == 'last':
            ax2 = ax.twinx()
            if j >= 0:
                ax2.set_ylabel(f'\nwMax\n= {j:.0f}', fontsize=12)
            else:
                ax2.set_ylabel('\nPoisson\nDispersion', fontsize=12)
            ax2.set_yticks([])


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with branch multiplier')
    parser.add_argument('input', type=str, help='Summary file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', nargs='+', type=float,
        help='wMax values to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    plot_wmax_pval(args.input, args.output, args.wMax)