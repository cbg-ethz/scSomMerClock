#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *


def generate_pval_plot_clock(args):
    df = pd.read_csv(args.input, sep='\t', index_col=0)
    df.drop([-1], inplace=True)
    df.drop(['R^2_pVal'], axis=1, inplace=True)

    cell_no = int(re.search('_bulk(\d+)_', args.input).group(1))

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))


    plot_neutralitytest(df['area_pVal'], ax)
    plot_mobster(df['s_Bayes'], ax2)

    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.75)
    fig2.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.75)
    if args.output:
        fig.savefig(args.output + '_neutralitytest.png', dpi=300)
        fig2.savefig(args.output + '_mobster.png', dpi=300)
    else:
        plt.show()
    plt.close()


def plot_neutralitytest(S, ax):
    color = '#B2DF8A'

    data = S.values.astype(float)
    dp = sns.histplot(data,
        element='bars', stat='probability', kde=False,
        binwidth=0.05, binrange=(0, 1),
        kde_kws={'cut': 0, 'clip': (0, 1)},
        line_kws={'lw': 3}, legend=False, ax=ax,
        color=color,
    )
    ax.set_xlim((0, 1))

    # Add rugplots
    add_rugs(data, offset=0, ax=ax, color=color)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('P-value')
    ax.set_ylabel('Probability')

    ax2 = ax.twinx()
    ax2.set_ylabel('\nneutralitytest', fontsize=12)
    ax2.set_yticks([])


def plot_mobster(S, ax):
    color = '#33A02C'

    data = S.values.astype(float)
    dp = sns.histplot(data,
        element='bars', stat='probability', kde=False, fill=True,
        kde_kws={'cut': 0}, line_kws={'lw': 3}, legend=False, ax=ax,
        color=color, log_scale=(False, False)
    )
    ax.set_ylim((0, 0.5))

    add_rugs(data, offset=0, ax=ax, color=color)
    ax.axvline(np.mean(data[~np.isnan(data)]), ls='--', color='black', lw=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Selective advantage s')
    ax.set_ylabel('Probability')

    ax2 = ax.twinx()
    ax2.set_ylabel('\nMobster', fontsize=12)
    ax2.set_yticks([])


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with branch multiplier')
    parser.add_argument('input', type=str, help='Summary file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_clock(args)