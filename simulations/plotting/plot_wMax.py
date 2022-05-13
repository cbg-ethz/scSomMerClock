#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
from matplotlib.lines import Line2D
import pandas as pd

from defaults import *


def generate_wmax_plot(args):
    df = pd.DataFrame(columns=['ADO rate', 'wMax', 'tree', 'significant'])
    cmap =cm.get_cmap('viridis_r')

    for res_file in os.listdir(args.input):
        if  not res_file.startswith('res_clock0'):
            continue

        if args.bulk and not 'bulk' in res_file:
            continue
        elif not args.bulk and 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))

        new_df = pd.read_csv(os.path.join(args.input, res_file), sep='\t', index_col=0)
        avg_row = new_df.loc[-1]
        new_df.drop(-1, inplace=True)

        for col, content in new_df.items():
            if not 'poissonTree' in col or not 'p-value' in col:
                continue
            try:
                wMax = int(re.search('_wMax([\.\d]*)\.', col).group(1))
            except AttributeError:
                continue
            tree = col.split('.')[-1]
            if not tree in METHODS:
                tree = col.split('_')[-1]

            if not tree in args.method:
                continue
            sig = (content < 0.05).mean()

            df.loc[df.shape[0]] = [ADO, wMax, tree, sig]
            # if ADO * 1000 == wMax:
            #     df.loc[df.shape[0]] = [ADO, -1, tree, sig]

        if (df[(df['ADO rate'] == ADO) & (df['wMax'] == 1)].size == 0) \
                & ((df['ADO rate'] == ADO).sum() > 0):
            for tree in ['cellcoal', 'cellphy', 'scite']:
                df.loc[df.shape[0]] = [ADO, 1, tree, 1]

    single_plot = len(args.method) == 1
    trees = df['tree'].unique()

    df['ADO rate'] /= 2
    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig, axes = plt.subplots(nrows=trees.size, ncols=1,
            figsize=(6, 2 * trees.size))
    axes = np.reshape(axes, (trees.size, 1))

    for i, tree in enumerate(trees):
        ax = axes[i][0]
        sns.lineplot(data=df[df['tree'] == tree], x='ADO rate', y='significant',
            hue='wMax', ax=ax, palette=cmap, marker='o', legend='full'
        )
        ax.set_ylabel(r'P-values $\leq$ 0.05 [%]')
        ax.axhline(0.05, ls=(0, (1, 4)), color='red', lw=2)
        ax.set_ylim((-0.05, 1.05))
        ax.set_xlabel('FN rate')
        ax.text(0.0125 / 2, 0.07, '0.05', color='red', ha='right', va='center')

        if i != trees.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if not single_plot:
            ax2 = ax.twinx()
            ax2.set_ylabel(METHODS[tree])
            ax2.set_yticks([])

        ax.get_legend().remove()

    if args.legend:
        legend_data = ax.get_legend_handles_labels()
        generate_legend_plot(*legend_data, args.output)

    fig.tight_layout()
    if args.output:
        fig.savefig(args.output, dpi=DPI)
    else:
        plt.show()
    plt.close()


def generate_legend_plot(handles, labels, output, circle=True):
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.grid(False)
    ax.axis('off')
    final_handles = []
    for handle in handles:
        if circle:
            final_handles.append(Line2D(
                [0], [0], marker='o', color=handle.get_color(), markersize=8, lw=0))
        else:
            handle.set_lw(2)
            final_handles.append(handle)
    ax.legend(final_handles, labels, ncol=2, frameon=True, title=r'$\bf{w_{max}}$')

    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    if output:
        fig.savefig(os.path.splitext(output)[0] + '_legend.png', dpi=DPI)
    else:
        plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, default='.',
        help='Base dir for results.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-b', '--bulk', action='store_true',
        help='Consider only bulk file.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        default = ['cellcoal', 'cellphy', 'scite'],
        choices=['cellcoal', 'cellphy', 'scite'],
        help='Which tree inference method to plot.')
    parser.add_argument('-l', '--legend', action='store_true',
        help='Plot legend as separate figure.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_wmax_plot(args)