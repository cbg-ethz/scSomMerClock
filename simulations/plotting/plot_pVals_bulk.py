#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *


def generate_pval_plot_bulk(args):
    df = pd.DataFrame(
        columns=['area_pVal', 's_Bayes', 'clones_Bayes', 'aff. cells', 'Amplifier'])
    for res_file in sorted(os.listdir(args.input)):
        if not res_file.startswith('res_clock') or not 'bulk' in res_file:
            continue
        bulk_file = os.path.join(args.input, res_file)
        df_new = pd.read_csv(bulk_file, sep='\t', index_col='run')

        df_new.drop([-1], inplace=True)
        df_new.drop(['R^2_pVal'], axis=1, inplace=True)

        if res_file.startswith('res_clock00'):
            ampl = 1
        else:
            ampl = float(re.search('res_clock(\d[\.\d]*)x', res_file).group(1))
        if args.amplifier and ampl not in args.amplifier:
            continue
        df_new['Amplifier'] = ampl

        print(res_file)
        for cl, data in df_new.groupby('clones_Bayes'):
            if ampl > 1:
                data = data[data['aff. cells'] > 30]
            print(f'\tC{cl:.0f}: n={data.shape[0]: >3}\t{data["s_Bayes"].mean():.4f}')

        df = df.append(df_new, ignore_index=True)

    cell_no = int(re.search('_bulk(\d+)_', bulk_file).group(1))

    ampl_vals = df['Amplifier'].unique()
    ampl_vals.sort()

    clone_sizes = np.array(args.clone_size)

    single_plot = clone_sizes.size < 3 and ampl_vals.size == 1

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
        fig2, axes2 = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
        axes = np.reshape(axes, (1, 1))
        axes2 = np.reshape(axes2, (1, 1))
    elif clone_sizes.size == 2:
        fig, axes = plt.subplots(nrows=ampl_vals.size, ncols=1,
            figsize=(3, ampl_vals.size + 2))
        fig2, axes2 = plt.subplots(nrows=ampl_vals.size, ncols=1,
            figsize=(3, ampl_vals.size + 2))
        axes = np.reshape(axes, (ampl_vals.size, 1))
        axes2 = np.reshape(axes2, (ampl_vals.size, 1))
    else:
        fig, axes = plt.subplots(nrows=ampl_vals.size, ncols=clone_sizes.size,
            figsize=(3 * clone_sizes.size, ampl_vals.size + 2))
        fig2, axes2 = plt.subplots(nrows=ampl_vals.size, ncols=clone_sizes.size,
            figsize=(3 * clone_sizes.size, ampl_vals.size + 2))
        axes = np.reshape(axes, (ampl_vals.size, clone_sizes.size))
        axes2 = np.reshape(axes2, (ampl_vals.size, clone_sizes.size))

    for i, clone_size_max in enumerate(clone_sizes):
        if clone_size_max == 0:
            df_plot = df
        elif clone_sizes.size == 2:
            min_cells = clone_sizes[0] * args.total_cells / 100
            max_cells = clone_sizes[1] * args.total_cells / 100
            df_plot = df[((df['aff. cells'] >= min_cells) \
                    & (df['aff. cells'] < max_cells)) \
                | (df['Amplifier'] == 1)]
        else:
            clone_size_min = clone_sizes[i - 1]
            min_cells = clone_size_min * args.total_cells / 100
            max_cells = clone_size_max * args.total_cells / 100
            df_plot = df[((df['aff. cells'] >= min_cells) \
                    & (df['aff. cells'] < max_cells)) \
                | (df['Amplifier'] == 1)]

        if df_plot.size == 0:
            print('!WARNING - No run with amplified clone sizes: '\
                f'[{clone_size_min}, {clone_size_max}]!')
            continue

        if clone_sizes.size == 2:
            fig3 = plot_mobster_box(df_plot)
            plot_neutralitytest(df_plot, axes[:, i], (i, 1))
            plot_mobster(df_plot, axes2[:, i], (i, 1))
            break
        else:
            plot_neutralitytest(df_plot, axes[:, i], (i, clone_sizes.size))
            plot_mobster(df_plot, axes2[:, i], (i, clone_sizes.size))

        if not single_plot:
            if clone_size_max == 0:
                header = 'All data'
            else:
                header = 'Amplified clone size:\n' \
                    f'({clone_size_min}, {clone_size_max}] %'
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
        fig.subplots_adjust(left=0.2, right=0.85, bottom=0.15, top=0.95,
            hspace=0.5, wspace=0.5)
        fig2.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.85,
            hspace=0.5, wspace=0.5)

    if args.output:
        fig.savefig(args.output + '_neutralitytest.png', dpi=300)
        fig2.savefig(args.output + '_mobster.png', dpi=300)
        if clone_sizes.size == 2:
            fig3.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.85)
            fig3.savefig(args.output + '_mobster_box.png', dpi=300)
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
        ax.set_ylim((0, 0.75))

        if ampl == 1:
            ax.set_yticks([0.05, 0.25, 0.5])
        else:
            ax.set_yticks([0.25, 0.5])

        add_rugs(data, offset=0, ax=ax, color=colors['neutrality'], height=0.05)

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
            if ampl == 1:
                ax2.set_ylabel(f'\nClock')
            else:
                ax2.set_ylabel(f'\nAmplifier:\n{ampl:.0f}')
            ax2.set_yticks([])


def plot_mobster_box(df):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    ax2 = ax.twiny()

    bar_vals = []
    for ampl, ampl_df in df.groupby('Amplifier'):
        for cl, cl_df in ampl_df.groupby('clones_Bayes'):
            bar_vals.append([ampl, cl, cl_df.shape[0] / ampl_df.shape[0] * 100])
    bar_df = pd.DataFrame(bar_vals, columns=['Amplifier', 'clones_Bayes', 'sum'])

    sns.barplot(data=bar_df, y="Amplifier", x="sum", hue='clones_Bayes',
        hue_order=[0, 1, 2], palette=colors, ax=ax2, dodge=True, alpha=.2,
        orient='h')
    sns.stripplot(data=df, y="Amplifier", x="s_Bayes", hue='clones_Bayes',
        hue_order=[0, 1, 2], palette=colors, ax=ax, dodge=True, linewidth=1,
        jitter=0.25, alpha=.4, orient='h', size=2)
    sns.boxplot(data=df, y="Amplifier", x="s_Bayes", hue='clones_Bayes',
        hue_order=[0, 1, 2], palette=colors, ax=ax, fliersize=2, orient='h',
        linewidth=1)

    ax2.set_xticks(np.array([100/6]) * np.arange(-2, 7, 1))
    ax2.set_xticklabels(['', '', 0, 16.7, 33.3, 50, 66.7, 83.3, 100])
    ax2.grid(False)
    ax.set_xlim((-5, 15))

    y_labels = []
    for i in sorted(df['Amplifier'].unique()):
        if i == 1:
            y_labels.append('Clock')
        else:
            y_labels.append(f'{i: >2.0f}')

    ax.set_ylabel('Amplifier')
    ax.set_yticklabels(y_labels)
    ax.set_xlabel('Fitness coefficient s')
    ax2.set_xlabel('Frequency')

    handles, labels = ax.get_legend_handles_labels()
    ax2.get_legend().remove()
    ax.legend(handles[:3], labels[:3], title='Clones', facecolor='white',
        framealpha=1)

    return fig


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

        data = df[df['Amplifier'] == ampl][['s_Bayes', 'clones_Bayes']]
        runs = data.shape[0]
        data = data.dropna()
        # data = data[(data['s_Bayes'] > -10) & (data['s_Bayes'] < 100)]

        dp = sns.histplot(data, x='s_Bayes', hue='clones_Bayes',
            element='bars', stat='probability',
            binwidth=0.1,
            kde=False, fill=True,
            palette=colors, log_scale=(False, False), legend=False, ax=ax,
        )

        ax.set_xlim((-0.5, 15))
        # ax.set_ylim((0, 0.7))
        ax.set_yticks([0.25, 0.5])

        for cl, cl_data in data.groupby('clones_Bayes'):
            add_rugs(cl_data['s_Bayes'], offset=0, ax=ax, color=colors[cl])
            if cl > 0:
                ax.axvline(cl_data['s_Bayes'].mean(), ls='--', color=colors[cl],
                    alpha=0.9, lw=1)

        # ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if col == 0:
            ann_text = f'n = {data.shape[0]} / {runs}'
        else:
            ann_text = f'n = {data.shape[0]}'
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
    parser.add_argument('-t', '--total_cells', type=int, default=100,
        help='Number of simulated cells. Default = 100')
    parser.add_argument('-c', '--clone_size', nargs='+',
        default = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        #default = [0, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95],
        type=float, help='Amplified clone size subsets. Default = [0, 10, 20, 30, 40, 50].')
    parser.add_argument('-a', '--amplifier', nargs='+', type=float,
        default=[1, 2, 5, 10, 20],
        help='Amplifier values to plot. Clock = 1. Default = [1, 5, 10, 20].')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_bulk(args)