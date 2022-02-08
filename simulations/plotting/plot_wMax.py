#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd

from defaults import *


def plot_wmax_pval2(in_dir, out_file=''):
    df = pd.DataFrame(columns=['ADO rate', 'wMax', 'False pos. [%]'])
    for res_dir in os.listdir(in_dir):
        if not res_dir.startswith('res_clock0') or 'bulk' in res_dir:
            continue
        ADO = float(re.search('WGA(0[\.\d]*)', res_dir).group(1))
        filter_dir = os.path.join(in_dir, res_dir, 'minDP5-minGQ1')
        for pTree_dir in os.listdir(filter_dir):
            if not pTree_dir.startswith('poissonTree'):
                continue
            wMax = int(re.search('wMax(\d+)', pTree_dir).group(1))
            wMax_file = os.path.join(filter_dir, pTree_dir, 'poissonTree.summary.tsv')
            with open(wMax_file, 'r') as f:
                res_raw = f.read().strip().split('\n')
            sig_raw =  res_raw[-1].split('\t')[-1].split('/')
            sig = 1 - float(sig_raw[0]) / float(sig_raw[1])

            if np.isclose(wMax, 1000 * ADO, atol=5):
                df.loc[df.shape[0]] = [ADO, r'$(\gamma + \beta)$ 1000', sig]
            else:
                df.loc[df.shape[0]] = [ADO, wMax, sig]

    fig, ax = plt.subplots(figsize=(12, 9))
    sns.lineplot(data=df, x='ADO rate', y='False pos. [%]', hue='wMax',
        palette={101: '#FF7F00', r'$(\gamma + \beta)$ 1000': '#377DB8',
            999: '#E41A1A'}, ax=ax, lw=2
    )

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_wmax_pval(in_dir, out_file='', bulk=False):
    df = pd.DataFrame(columns=['ADO rate', 'wMax', 'tree', 'significant'])
    cmap =cm.get_cmap('viridis_r')

    for res_file in os.listdir(in_dir):
        if  not res_file.startswith('res_clock0'):
            continue

        if bulk and not 'bulk' in res_file:
            continue
        elif not bulk and 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))

        new_df = pd.read_csv(os.path.join(in_dir, res_file), sep='\t', index_col=0)
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
            if not tree in vis_names:
                tree = col.split('_')[-1]
            sig = (content < 0.05).mean()

            df.loc[df.shape[0]] = [ADO, wMax, tree, sig]
            # if ADO * 1000 == wMax:
            #     df.loc[df.shape[0]] = [ADO, -1, tree, sig]

        if (df[(df['ADO rate'] == ADO) & (df['wMax'] == 1)].size == 0) \
                & ((df['ADO rate'] == ADO).sum() > 0):
            for tree in ['cellcoal', 'cellphy', 'scite']:
                df.loc[df.shape[0]] = [ADO, 1, tree, 1]

    trees = df['tree'].unique()
    fig, axes = plt.subplots(nrows=trees.size, ncols=1,
        figsize=(6, 2 * trees.size))

    for i, tree in enumerate(trees):
        ax = axes[i]
        if i == trees.size - 1:
            legend_type = 'full'
        else:
            legend_type = False

        sns.lineplot(data=df[df['tree'] == tree], x='ADO rate', y='significant',
            hue='wMax', ax=ax, lw=2, palette=cmap, marker='o',
            legend=legend_type
        )
        ax.set_ylabel(r'p-values $\leq$ 0.05 [%]')
        ax.axhline(0.05, ls='--', color='red', lw=2)
        ax.set_ylim((-0.05, 1.05))

        if i != trees.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        ax2 = ax.twinx()
        ax2.set_ylabel(vis_names[tree])
        ax2.set_yticks([])

    handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()
    ax.legend(handles, labels, ncol=2, bbox_to_anchor=(1.1, 1), # right, top
        frameon=True, title=r'$w_{max}$')

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.65, top=0.95, wspace=0.75)
    if out_file:
        fig.savefig(out_file, dpi=DPI)
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    plot_wmax_pval(args.input, args.output, args.bulk)