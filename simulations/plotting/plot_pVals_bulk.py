#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *


def generate_pval_plot_bulk(args):
    df = pd.DataFrame(columns=['area_pVal', 's_Bayes', 'aff. cells', 'Amplifier'])
    for res_file in os.listdir(args.input):
        if not res_file.startswith('res_clock') or not 'bulk' in res_file:
            continue
        bulk_file = os.path.join(args.input, res_file)
        df_new = pd.read_csv(bulk_file, sep='\t', index_col='run')
        df_new.drop([-1], inplace=True)
        df_new.drop(['R^2_pVal'], axis=1, inplace=True)

        if res_file.startswith('res_clock0'):
            ampl = 1
        else:
            ampl = float(re.search('res_clock(\d[\.\d]*)x', res_file).group(1))
        if args.amplifier and ampl not in args.amplifier:
            continue
        df_new['Amplifier'] = ampl

        df = df.append(df_new, ignore_index=True)
    cell_no = int(re.search('_bulk(\d+)_', bulk_file).group(1))

    ampl_vals = df['Amplifier'].unique()
    ampl_vals.sort()

    min_cells = np.array(args.min_cell)

    single_plot = min_cells.size == 1 and ampl_vals.size == 1

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
        fig2, axes2 = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig, axes = plt.subplots(nrows=ampl_vals.size, ncols=min_cells.size,
            figsize=(3 * min_cells.size, ampl_vals.size + 2))
        fig2, axes2 = plt.subplots(nrows=ampl_vals.size, ncols=min_cells.size,
            figsize=(3 * min_cells.size, ampl_vals.size + 2))
    axes = np.reshape(axes, (ampl_vals.size, min_cells.size))
    axes2 = np.reshape(axes2, (ampl_vals.size, min_cells.size))

    for i, min_cell in enumerate(min_cells):
        df_plot = df[(df['aff. cells'] >= min_cell) | (df['Amplifier'] == 1)]
        if df_plot.size == 0:
            print(f'!WARNING - No run with {min_cell} cells affected!')
            continue

        plot_neutralitytest(df_plot, axes[:, i], (i, min_cells.size))
        plot_mobster(df_plot, axes2[:, i], (i, min_cells.size))

        if not single_plot:
            header = r'$\geq$' + f'{min_cell}/{cell_no} cells'
            axes[0,i].annotate(header, xy=(0.5, 1.15), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
            axes2[0,i].annotate(header, xy=(0.5, 1.15), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

    if single_plot:
        fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.75)
        fig2.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.75)
    else:
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.85,
            hspace=0.5, wspace=0.5)
        fig2.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.85,
            hspace=0.5, wspace=0.5)

    if args.output:
        fig.savefig(args.output + '_neutralitytest.png', dpi=300)
        fig2.savefig(args.output + '_mobster.png', dpi=300)
    else:
        plt.show()
    plt.close()


def plot_neutralitytest(df, axes, col):
    ampl_vals = sorted(df['Amplifier'].unique())

    for i, ampl in enumerate(ampl_vals):
        ax = axes[i]

        if ampl == 1 and col[0] != 0:
            ax.grid(False)
            ax.axis('off')
            if col[0] == col[1] - 1:
                ax2 = ax.twinx()
                ax2.grid(False)
                for spine in ax.spines:
                    ax2.spines[spine].set_visible(False)
                ax2.set_xticks([])
                ax2.set_yticks([])
                ax2.set_ylabel(f'\nClock')
            continue

        data = df[df['Amplifier'] == ampl]['area_pVal'].values.astype(float)
        dp = sns.histplot(data,
            element='bars', stat='probability', kde=False, binwidth=0.05,
            binrange=(0, 1), color=colors['neutrality'], legend=False, ax=ax,
        )
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 0.6))

        if ampl == 1:
            ax.set_yticks([0.05, 0.25, 0.5])
        else:
            ax.set_yticks([0.25, 0.5])

        add_rugs(data, offset=0, ax=ax, color=colors['neutrality'])

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.annotate(f'n = {data.size:.0f}', xy=(0.95, 0.75), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='top', bbox=bbox_props)

        if col[0] == np.floor((col[1] - 1) / 2) and i == len(ampl_vals) - 1:
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')

        if col[0] == 0 and i == np.floor((len(ampl_vals) - 1) / 2):
            ax.set_ylabel('Probability')
        else:
            ax.set_ylabel('')

        if col[0] == col[1] - 1:
            ax2 = ax.twinx()
            ax2.set_ylabel(f'\nAmplifier:\n{ampl:.0f}x')
            ax2.set_yticks([])


def plot_mobster(df, axes, col):
    ampl_vals = sorted(df['Amplifier'].unique())

    for i, ampl in enumerate(ampl_vals):
        ax = axes[i]

        if ampl == 1 and col[0] != 0:
            ax.grid(False)
            ax.axis('off')
            if col[0] == col[1] - 1:
                ax2 = ax.twinx()
                ax2.grid(False)
                for spine in ax.spines:
                    ax2.spines[spine].set_visible(False)
                ax2.set_xticks([])
                ax2.set_yticks([])
                ax2.set_ylabel(f'\nClock')
            continue

        data = df[df['Amplifier'] == ampl]['s_Bayes'].values.astype(float)
        data = np.clip(data, -10, 100)

        dp = sns.histplot(data,
            element='bars', stat='probability', kde=False, fill=True,
            color=colors['mobster'], log_scale=(False, False), legend=False, ax=ax,
        )
        ax.set_xlim((-2, 20))
        ax.set_ylim((0, 0.7))
        ax.set_yticks([0.25, 0.5])

        add_rugs(data, offset=0, ax=ax, color=colors['mobster'])
        ax.axvline(np.mean(data[~np.isnan(data)]), ls='--', color='black',
            alpha=0.75, lw=1)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if col == 0:
            ann_text = f'n = {np.isfinite(data).sum()} / {data.size}'
        else:
            ann_text = f'n = {np.isfinite(data).sum()}'
        ax.annotate(ann_text, xy=(0.95, 0.75), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='top', bbox=bbox_props)

        if col[0] == np.floor((col[1] - 1) / 2) and i == len(ampl_vals) - 1:
            ax.set_xlabel('Selective advantage s')
        else:
            ax.set_xlabel('')

        if col[0] == 0 and i == np.floor((len(ampl_vals) - 1) / 2):
            ax.set_ylabel('Probability')
        else:
            ax.set_ylabel('')

        if col[0] == col[1] - 1:
            ax2 = ax.twinx()
            ax2.set_ylabel(f'\nAmplifier:\n{ampl:.0f}x')
            ax2.set_yticks([])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Directory with summary files.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-c', '--min_cell', nargs='+', default = [1, 10, 20, 30, 40, 50],
        type=float, help='Min. #cells affected. Default = [1, 10, 20, 30, 40, 50].')
    parser.add_argument('-a', '--amplifier', nargs='+', default=[], type=float,
        help='Clock and amplifier values to plot. Clock = 1. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_bulk(args)