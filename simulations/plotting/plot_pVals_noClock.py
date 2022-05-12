#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *
CELL_NO = 30


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
            dpi=DPI)
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
                if tree not in args.method:
                    continue
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
    cols =  ['ADO', 'amplifier', 'aff. cells', 'method', 'P-value']
    df = pd.DataFrame(columns=cols)

    for res_file in sorted(os.listdir(args.input)):
        if not res_file.startswith('res_clock') or 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO not in args.ADO:
            continue

        ampl = float(re.search('res_clock(\d[\.\d]*)', res_file).group(1))
        if ampl == 0:
            ampl = 1

        if ampl not in args.amplifier:
            continue

        print(f'Including summary file: {res_file}')

        df_new = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_new.drop([-1], inplace=True)

        if ampl > 1:
            df_new = df_new[
                (df_new['aff. cells'] >= args.total_cells * args.clone_size[0]) \
                & (df_new['aff. cells'] <= args.total_cells * args.clone_size[1])]

        for col in df_new.columns:
            if col.startswith('p-value'):
                if col.endswith('poissonDisp'):
                    tree = 'poissonDisp'
                    wMax = args.wMax
                elif 'paup_dir' in col:
                    tree = 'PAUP*'
                    wMax = args.wMax
                else:
                    tree = re.split('[\._]', col)[-1]
                    if ampl == 1:
                        try:
                            wMax = int(re.search('wMax(\d+)', col).group(1))
                        except AttributeError:
                            wMax = int(re.split('[\._]', col)[2])
                            if wMax == 1:
                                wMax = args.wMax
                    else:
                        wMax = int(re.search('wMax(\d+)', col).group(1))

                    if wMax != args.wMax:
                        continue

                if tree not in args.method:
                    continue

                if ampl == 1:
                    df_subset = df_new[[col]]
                else:
                    df_subset = df_new[['aff. cells', col]]

                df_subset.insert(0, 'ADO', ADO)
                df_subset.insert(0, 'method', tree)
                df_subset.insert(0, 'amplifier', ampl)
                df_subset = df_subset.rename({col: 'P-value'}, axis=1)

                df = df.append(df_subset, ignore_index=True)

    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()

    ampl_vals = df['amplifier'].unique()
    ampl_vals.sort()

    col_no = ADO_vals.size
    row_no = ampl_vals.size

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no + 1, row_no + 1))
    axes = np.reshape(axes, (row_no, col_no))

    for i, ampl_val in enumerate(ampl_vals):
        df_plot = df[df['amplifier'] == ampl_val]
        plot_pVal_dist(df_plot, ADO_vals, axes[i], (i, row_no))
        if ampl_val == 1:
            row_title = 'Clock'
        else:
            row_title = f'Amplifier:\n{ampl_val:.0f}x'
        add_row_header(axes[i, -1], row_title)

    if col_no > 1:
        for j, ADO_val in enumerate(ADO_vals):
            add_col_header(axes[0, j], f'FN rate: {ADO_val / 2}')

    fig.tight_layout()
    if args.output:
        if not args.output.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            args.output += '.png'
        fig.savefig(args.output, dpi=DPI)
    else:
        plt.show()
    plt.close()


def plot_pVal_dist(df, ADO_vals, axes, row):
    for i, ADO_val in enumerate(ADO_vals):
        ax = axes[i]
        data = df[df['ADO'] == ADO_val]

        dp = sns.histplot(data, x='P-value', hue='method', ax=ax, **HIST_DEFAULT)
        # # import pdb; pdb.set_trace()
        # dp = sns.violinplot(x='wMax', y='P-value', data=df, hue='method',
        #     cut=0, split=True, inner='point',
        #     palette=colors,
        #     hue_order=hue_order,
        #     ax=ax,
        # )

        ax.set_xlim((0, 1.01))
        ax.set_xticks([0.0, 0.5, 1])
        ax.set_ylim((0, 1))
        ax.set_yticks([0.0, 0.5, 1])

        # Add rugplots
        k = 0
        for method in HUE_ORDER:
            rug_data = data[data['method'] == method]['P-value'].values
            if rug_data.size == 0:
                continue
            add_rugs(rug_data, offset=k, ax=ax, color=COLORS[method])
            k += 1

        ax.annotate(f'n = {data["method"].value_counts().min():.0f}',
            xy=(0.9, 0.825), xytext=(0, 5),  xycoords='axes fraction',
            textcoords='offset points', ha='right', va='top', bbox=bbox_props)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #
        if row[0] < row[1] - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        # First col
        if i == 0:
            if row[0] != np.floor(row[1] / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel('Probability')
        else:
            ax.set_ylabel(None)
            ax.set_yticklabels([])

        # Middle col
        if i == np.floor(ADO_vals.size / 2) and row[0] == row[1] - 1:
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')


def parse_args():
    parser = argparse.ArgumentParser(
        'Generate p-value plots for CellCoal simulations with branch multiplier')
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', type=float, default=400,
        help='wMax value to plot. Default = 400.')
    parser.add_argument('-a', '--ADO', nargs='+', type=float,
        default=[0, 0.2, 0.4], help='ADO values to plot. Default = [0, 0.2, 0.4].')
    parser.add_argument('-amp', '--amplifier', nargs='+', default=[2, 5, 10],
        type=float, help='Amplifier value to plot. Default = [2, 5, 10]')
    parser.add_argument('-c', '--clone_size', default = [0.1, 0.9],
        type=float, help='Amplified clone size subsets. Default = [0.1, 0.9].')
    parser.add_argument('-t', '--total_cells', type=int, default=30,
        help='Total number of simualted cells. Default = 30.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy', 'scite', 'poissonDisp'],
        default=['cellcoal', 'cellphy', 'poissonDisp'],
        help='Method to plot. Default = all.')
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