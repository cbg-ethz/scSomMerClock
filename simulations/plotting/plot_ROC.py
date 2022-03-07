#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

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


def generate_ROC(args):
    cols =  ['aff. cells', 'method', 'wMax', 'predicted']
    df = pd.DataFrame(columns=cols)
    for res_file in os.listdir(args.input):
        if not res_file.startswith('res_clock') or 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO != args.ADO:
            continue

        if res_file.startswith('res_clock0'):
            clock = True
        else:
            ampl = float(re.search('res_clock(\d[\.\d]*)x', res_file).group(1))
            if not ampl == args.amplifier:
                continue
            clock = False

        print(f'Including summary file: {res_file}')

        df_summary = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_summary.drop([-1], inplace=True)
        df_summary = df_summary.T.drop_duplicates().T

        rel_cols = [i for i in df_summary.columns if i.startswith('p-value')]
        for rel_col in rel_cols:
            if clock:
                new_df = pd.DataFrame(
                    np.where(df_summary[rel_col] < 0.05, 'FP', 'TN'),
                    columns=['predicted'])
                new_df['aff. cells'] = np.inf
            else:
                new_df = pd.DataFrame(
                    np.where(df_summary[rel_col] < 0.05, 'TP', 'FN'),
                    columns=['predicted'])
                new_df['aff. cells'] = df_summary['aff. cells']

            if rel_col[8:].startswith('poissonTree'):
                method = rel_col.split('.')[-1]
                if method not in colors:
                    method = rel_col.split('_')[-1]
                new_df['wMax'] = int(re.search('wMax(\d+)', rel_col).group(1))
                new_df['method'] = method
            elif rel_col[8:].startswith('poissonDisp'):
                continue
                # new_df['wMax'] = -1
                # new_df['method'] = 'Poisson Dispersion'
            else:
                raise TypeError(f'Unknown column type: {rel_col}')

            df = df.append(new_df, ignore_index=True)

    assert df['predicted'].unique().size == 4

    min_cells = np.array(args.min_cell)
    single_plot = min_cells.size == 1

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig, axes = plt.subplots(nrows=1, ncols=min_cells.size,
            figsize=(3 * min_cells.size, 3))
    axes = np.reshape(axes, (min_cells.size))

    for i, min_cell in enumerate(min_cells):
        df_plot = df[df['aff. cells'] >= min_cell]
        if df_plot.size == 0:
            print(f'!WARNING - No run with {min_cell} cells affected!')
            continue

        plot_ROC(df_plot, axes[i], (i, min_cells.size))

        if not single_plot:
            header = r'$\geq$' + f'{min_cell} cells'
            axes[i].annotate(header, xy=(0.5, 1.1), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

    if single_plot:
        fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.75)
    else:
        fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.85,
            hspace=0.5, wspace=0.5)

    if args.output:
        fig.savefig(args.output, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_ROC(df_in, ax, col_id):
    for method, df in df_in.groupby('method'):
        if method != 'Poisson Dispersion':
            rec = []
            prec = []
            w_max = []
            for wMax, df_w in df.groupby('wMax'):
                counts = df_w['predicted'].value_counts()
                vals = df_w['predicted'].value_counts()
                if vals.size < 3:
                    continue
                if not 'FP' in vals:
                    vals.FP = 0

                recall = vals.TP / (vals.TP + vals.FN) # Recall/TRP
                FPR = vals.FP / (vals.FP + vals.TN)
                precision = vals.TP / (vals.TP + vals.FP)

                w_max.append(wMax)
                prec.append(precision)
                rec.append(recall)

        else:
            vals = df['predicted'].value_counts()
            rec = [vals.TP / (vals.TP + vals.FN)]
            prec = [vals.TP / (vals.TP + vals.FP)]

        sns.lineplot(x=rec, y=prec, ax=ax, marker='o', markersize=4,
                color=colors[method], alpha  = 0.75)

    ax.set_xlim((0, 1.))
    ax.set_ylim((0.5, 1.05))

    ax.set_ylabel('Precision')
    ax.set_xlabel('Recall')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Directory with summary files.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-c', '--min_cell', nargs='+', default = [1, 2, 3, 4, 5],
        type=float, help='Min. #cells affected. Default = [1, 2, 3, 4, 5].')
    parser.add_argument('-a', '--amplifier', default=2, type=float,
        help='Amplifier value to plot. Default = 2.')
    parser.add_argument('-do', '--ADO', default=0.4, type=float,
        help='Simulated ADO value to plot. Default = 0.4.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_ROC(args)