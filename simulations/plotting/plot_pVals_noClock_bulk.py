#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=1)


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


def generate_pval_plot_noClock(args):
    df = pd.DataFrame(columns=['area_pVal', 's_Bayes', 'aff. cells', 'Amplifier'])
    for res_file in os.listdir(args.input):
        if not res_file.startswith('res_clock') or not 'bulk' in res_file \
                or res_file.startswith('res_clock0'):
            continue
        bulk_file = os.path.join(args.input, res_file)
        df_new = pd.read_csv(bulk_file, sep='\t', index_col='run')
        df_new.drop([-1], inplace=True)
        df_new.drop(['R^2_pVal'], axis=1, inplace=True)

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
        data = df[df['Amplifier'] == ampl]['area_pVal'].values.astype(float)
        dp = sns.histplot(data,
            element='bars', stat='probability', kde=False, binwidth=0.05,
            binrange=(0, 1), color=colors['neutrality'], legend=False, ax=ax,
        )
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 0.6))

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
        data = df[df['Amplifier'] == ampl]['s_Bayes'].values.astype(float)
        data = np.clip(data, -10, 100)

        dp = sns.histplot(data,
            element='bars', stat='probability', kde=False, fill=True,
            color=colors['mobster'], log_scale=(False, False), legend=False, ax=ax,
        )
        ax.set_xlim((-2, 20))
        ax.set_ylim((0, 0.7))

        add_rugs(data, offset=0, ax=ax, color=colors['mobster'])
        ax.axvline(np.mean(data[~np.isnan(data)]), ls='--', color='black',
            alpha=0.75, lw=1)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.annotate(f'n = {data.size:.0f}', xy=(0.95, 0.75), xytext=(0, 5),
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
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with branch multiplier')
    parser.add_argument('input', type=str, help='Summary file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-c', '--min_cell', nargs='+', default = [1, 10, 20, 30, 40, 50],
        type=float, help='Min. #cells affected. Default = [1, 10, 20, 30, 40, 50].')
    parser.add_argument('-w', '--amplifier', nargs='+', default=[], type=float,
        help='wMax values to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_noClock(args)