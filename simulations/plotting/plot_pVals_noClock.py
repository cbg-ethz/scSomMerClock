#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection

from defaults import *


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
    ax.set_ylabel(f'Rate of samples [%]')
    ax.set_xlim((0, 30))
    ax.set_xlabel(f'#Cells affected by switch')
    ax.tick_params(axis='both', which='major')
    ax.tick_params(axis='both', which='minor')
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
    ax2.set_ylabel(f'1 - Cum. Sum(Rate of samples)')
    ax2.tick_params(axis='both', which='major')
    ax2.tick_params(axis='both', which='minor')
    ax2.spines['right'].set_color(col2)
    ax2.tick_params(axis='y', colors=col2)

    if out_file:
        fig.savefig(os.path.splitext(out_file)[0] + '_affectedCells.png',
            dpi=300)
    else:
        plt.show()
    plt.close()


def plot_sign_over_cells(df_in, axes, trees, col_no):
    cmap =cm.get_cmap('viridis_r')

    for i, tree in enumerate(trees):
        df_tree = df_in[df_in['Tree'] == tree]
        data = []
        for cells, cells_data in df_tree.groupby('aff. cells'):
            for wMax, wMax_data in cells_data.groupby('wMax'):
                if wMax == -1:
                    continue
                data.append([wMax, cells, (wMax_data['P-value'] < 0.05).mean()])
        df = pd.DataFrame(data, columns=['wMax', 'cells', 'significant'])

        if (df['wMax'] == 1).sum() == 0:
            df.loc[df.shape[0]] = [1, np.nan, np.nan]

        ax = axes[i]
        sns.lineplot(data=df, x='cells', y='significant', hue='wMax', ax=ax,
            palette=cmap, legend=False
        )
        ax.set_ylim((-0.01, 1.01))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < trees.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if col_no[0] == 0:
            if i != np.floor(trees.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel(r'P-values $\leq$ 0.05 [%] (Recall)')
        else:
            ax.set_ylabel('')

        if col_no[0] == np.floor(col_no[1] / 2) and i == trees.size - 1:
            ax.set_xlabel(f'# affected cells (n = {df["cells"].max() + 1:.0f})')
        else:
            ax.set_xlabel('')

        if col_no[0] == col_no[1] - 1:
            ax2 = ax.twinx()
            ax2.set_ylabel(vis_names[tree], fontsize=12)
            ax2.set_yticks([])


def generate_sign_over_cells(args):
    df = pd.DataFrame(columns=['wMax', 'Tree', 'aff. cells', 'P-value', 'Amplifier'])
    for res_file in os.listdir(args.input):
        if not res_file.startswith('res_clock') or 'bulk' in res_file \
                or res_file.startswith('res_clock0'):
            continue
        ampl = float(re.search('res_clock(\d[\.\d]*)x', res_file).group(1))
        df_summary = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_summary.drop([-1], inplace=True)
        df_summary = df_summary.T.drop_duplicates().T

        cols = [i for i in df_summary.columns if i.startswith('p-value')]
        for col in cols:
            if col[8:].startswith('poissonTree'):
                try:
                    wMax = int(re.search('wMax(\d+)', col).group(1))
                except AttributeError:
                    # No errors
                    wMax = args.wMax
                else:
                    if wMax not in args.wMax:
                        continue
                    wMax = [wMax]
                tree = col.split('.')[-1]
                if tree not in colors:
                    tree = col.split('_')[-1]
            elif col[8:].startswith('poissonDisp') or col[8:].startswith('paup'):
                continue
            else:
                raise TypeError(f'Unknown column type: {col}')

            df_new = df_summary.loc[:,[col, 'aff. cells']]
            df_new.columns = ['P-value', 'aff. cells']

            df_new['Amplifier'] = ampl
            df_new['Tree'] = tree
            for wMax_i in wMax:
                df_new['wMax'] = wMax_i
                df = df.append(df_new, ignore_index=True)

    ampl_vals = df['Amplifier'].unique()
    ampl_vals.sort()

    tree_vals = df['Tree'].unique()

    fig, axes = plt.subplots(nrows=tree_vals.size, ncols=ampl_vals.size,
        figsize=(3 * ampl_vals.size, tree_vals.size + 2))
    axes = np.reshape(axes, (tree_vals.size, ampl_vals.size))

    for i, ampl_val in enumerate(ampl_vals):
        df_plot = df[df['Amplifier'] == ampl_val]
        plot_sign_over_cells(df_plot, axes[:,i], tree_vals, (i, ampl_vals.size))

        axes[0, i].annotate(f'Amplifier: {ampl_val}x', xy=(0.5, 1.1),
            xytext=(0, 5), xycoords='axes fraction',
            textcoords='offset points', size='large', ha='center',
            va='baseline')

    fig.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9,
            hspace=0.25, wspace=0.25)

    if args.output:
        fig.savefig(args.output, dpi=DPI)
    else:
        plt.show()
    plt.close()


def generate_pval_plot_noClock(args):
    df_in = pd.read_csv(args.input, sep='\t', index_col=0)
    df_in.drop([-1], inplace=True)
    rel_cols = ['aff. cells'] + [i for i in df_in.columns if 'p-value' in i]
    df_in = df_in[rel_cols]

    # plot_affected_cells(df_in, args.output)

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

    if args.wMax:
        df = df[df['wMax'].isin(args.wMax)]

    min_cells = np.array(args.min_cell)
    wMax_vals = df['wMax'].unique()
    wMax_vals.sort()

    single_plot = wMax_vals.size == 1 and min_cells.size == 1

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig, axes = plt.subplots(nrows=wMax_vals.size, ncols=min_cells.size,
            figsize=(3 * min_cells.size, wMax_vals.size + 1))
    axes = np.reshape(axes, (wMax_vals.size, min_cells.size))

    for i, min_cell in enumerate(min_cells):
        if single_plot:
            col_type = 'only'
        elif i == 0:
            col_type = 'first'
        elif i == min_cells.size - 1:
            col_type = 'last'
        elif i == np.floor(min_cells.size / 2):
            col_type = 'middle'
        else:
            col_type = 'intermediate'

        df_plot = df[df['aff. cells'] >= min_cell]
        if df_plot.size == 0:
            print(f'!WARNING - No run with {min_cell} cells affected!')
            continue
        plot_pVal_dist(df_plot, wMax_vals, axes[:,i], col_type)

        if not single_plot:
            wMax_max = df_plot['wMax'].value_counts().index[0]
            n = ((df_plot['Tree'] == 'cellcoal') & (df_plot['wMax'] == wMax_max)) \
                .sum()
            header = r'$\geq$' + f'{min_cell} cells\n(n = {n})'
            axes[0,i].annotate(header, xy=(0.5, 1.1), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

    if single_plot:
        fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.75)
    else:
        fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
            hspace=0.5, wspace=0.5)

    if args.output:
        fig.savefig(args.output, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_pVal_dist(df, wMax, axes, column_type):
    hue_order = ['-', 'cellcoal', 'cellphy', 'scite']
    for i, j in enumerate(wMax):
        ax = axes[i]
        data = df[df['wMax'] == j]
        dp = sns.histplot(data, x='P-value', hue='Tree',
            element='poly', stat='probability', kde=False,
            common_norm=False, fill=False,
            binwidth=0.01, binrange=(0, 1), multiple='dodge',
            kde_kws={'cut': 0, 'clip': (0, 1)},
            line_kws={'lw': 3, 'alpha': 0.66},
            palette=colors, alpha=0.5,
            hue_order=['-', 'cellcoal', 'cellphy', 'scite'],
            legend=False, ax=ax,
        )

        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))

        # Add rugplots
        k = 0
        for method in hue_order:
            rug_data = data[data['Tree'] == method]['P-value'].values
            if rug_data.size == 0:
                continue

            segs = np.stack((np.c_[rug_data, rug_data],
                np.c_[np.zeros_like(rug_data) + 1 + RUG_HEIGHT*2 * k,
                        np.zeros_like(rug_data) + 1 + RUG_HEIGHT*2 * (k + 1)]),
                    axis=-1)
            lc = LineCollection(segs,
                transform=ax.get_xaxis_transform(),
                clip_on=False, color=colors[method], linewidth=0.05, alpha=0.75)
            ax.add_collection(lc)
            k += 1

        l = 0
        for tree, col in colors.items():
            df_sub = df[(df['Tree'] == tree) & (df['wMax'] == j)]
            if df_sub.size > 0:
                y_pval = (df_sub['P-value'] <= 0.05).mean()
                if y_pval > 0.33:
                    va = 'top'
                    y_dist = -0.02
                else:
                    va = 'bottom'
                    y_dist = 0.02

                ax.axhline(y_pval, ls='--', color=col, lw=1)
                # ax.axvline(0.05, ls='--', color='grey', lw=1)
                ax.text(0.05 + l * 0.15, y_pval + y_dist, f'{y_pval:.2f}',
                    color=col, ha='left', va=va, rotation=45)
                l += 1

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < wMax.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if column_type == 'first':
            if i != np.floor(wMax.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel('Probability')
        else:
            ax.set_ylabel('')

        if column_type == 'middle' and i == wMax.size - 1:
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')

        if column_type == 'last':
            ax2 = ax.twinx()
            if j >= 0:
                ax2.set_ylabel('\n' + r'$w_{max}=$ '+ f'\n{j:.0f}', fontsize=12)
            else:
                ax2.set_ylabel('\nPoisson\nDispersion', fontsize=12)
            ax2.set_yticks([])

        if column_type == 'only':
            ax.set_xlabel('P-value')
            ax.set_ylabel('Probability')


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with branch multiplier')
    parser.add_argument('input', type=str, help='Summary file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', nargs='+', type=float,
        default=[1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
        help='wMax values to plot. Default = all.')
    parser.add_argument('-c', '--min_cell', nargs='+',default = [1, 2, 3, 4, 5],
        type=float, help='Min. #cells affected. Default = [1, 2, 3, 4, 5].')
    parser.add_argument('-scd', '--sign_cell_dist', action='store_true',
        help='Plot sign. p-values per aff. cells dist.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.sign_cell_dist:
        generate_sign_over_cells(args)
    else:
        generate_pval_plot_noClock(args)