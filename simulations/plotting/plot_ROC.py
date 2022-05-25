#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
from matplotlib.colors import to_rgba

from defaults import *
CELL_NO = 30


def generate_PRC(args):
    cols =  ['ADO', 'aff. cells', 'method', 'wMax', 'predicted']
    df = pd.DataFrame(columns=cols)

    for res_file in sorted(os.listdir(args.input)):
        if not res_file.startswith('res_clock') or 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO not in args.ADO:
            continue

        if re.match('res_clock0+_', res_file):
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
                if method not in COLORS:
                    method = rel_col.split('_')[-1]
                if method not in args.method:
                    continue
                if ADO == 0 and clock:
                    try:
                        wMax = int(rel_col.split('_')[-1].split('.')[0])
                    except:
                        wMax = 1
                    if wMax == 0:
                        continue
                else:
                    wMax = int(re.search('wMax(\d+)', rel_col).group(1))
                new_df['wMax'] = wMax
                new_df['method'] = method
            elif rel_col[8:].startswith(('poissonDisp', 'paup')):
                continue
                # new_df['wMax'] = -1
                # new_df['method'] = 'Poisson Dispersion'
            else:
                raise TypeError(f'Unknown column type: {rel_col}')

            new_df['ADO'] = ADO
            df = df.append(new_df, ignore_index=True)

    assert df['predicted'].unique().size == 4

    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()
    col_no = ADO_vals.size

    clone_sizes = np.array(args.clone_size)
    if clone_sizes.size < 3:
        row_no = 1
    else:
        row_no = clone_sizes.size

    fig, axes = get_subplots(row_no, col_no)

    for i, clone_size_max in enumerate(clone_sizes):
        if clone_size_max == 0:
            df_plot = df
        elif clone_sizes.size == 2:
            min_cells = clone_sizes[0] * args.total_cells / 100
            max_cells = clone_sizes[1] * args.total_cells / 100
            df_plot = df[((df['aff. cells'] >= min_cells) \
                    & (df['aff. cells'] < max_cells)) |
                (df['aff. cells'] == np.inf)]
        else:
            clone_size_min = clone_sizes[i - 1]
            min_cells = clone_size_min * args.total_cells / 100
            max_cells = clone_size_max * args.total_cells / 100
            df_plot = df[((df['aff. cells'] >= min_cells) \
                    & (df['aff. cells'] < max_cells)) |
                (df['aff. cells'] == np.inf)]

        if df_plot.size == 0:
            print('!WARNING - No run with amplified clone sizes: '\
                f'[{clone_size_min}, {clone_size_max}]!')
            continue

        plot_curve(df_plot, axes[i], ADO_vals, (i, row_no), curve=args.curve)

        if clone_sizes.size == 2:
            break

    for i, ADO_val in enumerate(ADO_vals):
        if ADO_val == 0:
            col_title = 'No Errors'
        else:
            col_title = f'FN rate: {ADO_val / 2}'
        add_col_header(axes[0, i], col_title)

    plot_fig(fig, args.output)



def plot_curve(df_in, axes,ADO_vals, row, row_title='', curve='ROC'):
    wMax_map = {j: i for i, j in enumerate(sorted(df_in['wMax'].unique()))}
    cmap = cm.get_cmap('viridis_r', len(wMax_map))

    for j, ADO_val in enumerate(ADO_vals):
        df_ADO = df_in[df_in['ADO'] == ADO_val]
        ax = axes[j]

        for method, df in df_ADO.groupby('method'):
            if method != 'Poisson Dispersion':
                rec = [0, 1]
                prec = [1, 0]
                fprs = [0, 1]
                w_max = []
                for wMax, df_w in df.groupby('wMax'):
                    vals = df_w['predicted'].value_counts()
                    if vals.size < 3:
                        continue
                    if not 'FP' in vals:
                        vals.FP = 0
                    if not 'TN' in vals:
                        vals.TN = 0

                    recall = vals.TP / (vals.TP + vals.FN) # Recall/Sensitivity/TPR
                    fpr = vals.FP / (vals.FP + vals.TN) #
                    precision = vals.TP / (vals.TP + vals.FP)

                    w_max.append(wMax)
                    prec.append(precision)
                    rec.append(recall)
                    # <INFO> Max to avoid bug in plotting vertical line to (0|0)
                    fprs.append(max(fpr, 1e-4))
            else:
                vals = df['predicted'].value_counts()
                rec = [vals.TP / (vals.TP + vals.FN)]
                prec = [vals.TP / (vals.TP + vals.FP)]

            if curve == 'ROC':
                sns.lineplot(x=fprs, y=rec, ax=ax, color=COLORS[method],
                    markersize=0, alpha=0.75, lw=1)
            else:
                sns.lineplot(x=rec, y=prec, ax=ax, color=COLORS[method],
                    markersize=0, alpha=0.75, lw=1)

            for k, wMax in enumerate(w_max, 2):
                color = [i for i in to_rgba(COLORS[method])[:3]] + [0.75]
                if curve == 'ROC':
                    ax.plot(fprs[k], rec[k], marker='o', markersize=4,
                        mec=color, mfc=cmap(k - 2))
                    ax.plot([0, 1], [0, 1], color='grey', ls='--', lw=1)
                else:
                    ax.plot(rec[k], prec[k], marker='o', markersize=4,
                        mec=color, mfc=cmap(k - 2))


        format_ax(ax, row, (j, ADO_vals.size), curve == 'PRC', row_title)


def format_ax(ax, row, col, PRC, row_title):
    if PRC:
        ax.set_xlim((0.05, 1.05))
        ax.set_xticks([0, 0.5, 1])
        ax.set_ylim((0.4, 1.05))
        ax.set_yticks([0.5, 1])
    else:
        ax.set_xlim((-0.05, 1.05))
        ax.set_xticks([0, 0.5, 1])
        ax.set_ylim((0, 1.05))
        ax.set_yticks([0, 0.5, 1])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # not bottom row
    if row[0] < row[1] - 1:
        ax.set_xticklabels([])
        ax.set_xlabel(None)

    # First col
    if col[0] == 0:
        # middle row
        if row[0] == np.floor(row[1] / 2):
            if PRC:
                ax.set_ylabel('Precision')
            else:
                ax.set_ylabel('True positive rate')
        else:
            ax.set_ylabel(None)
    else:
        ax.set_ylabel(None)
        ax.set_yticklabels([])

    # Middle col, last row
    if col[0] == np.floor((col[1] - 1) / 2) and row[0] == row[1] - 1:
        if PRC:
            ax.set_xlabel('Recall')
        else:
            ax.set_xlabel('False positive rate')
    else:
        ax.set_xlabel('')

    # Last column
    if col[0] == col[1] - 1:
        add_row_header(ax, f'\n{row_title}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Directory with summary files.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-cu', '--curve', type=str, default='ROC',
        choices=['ROC', 'PRC'], help='Which curve to plot. Default = ROC.')
    parser.add_argument('-t', '--total_cells', type=int, default=30,
        help='Number of simulated cells. Default = 30')
    parser.add_argument('-c', '--clone_size', nargs='+',
        default = [10, 90],
        #default = [0, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95],
        type=float, help='Amplified clone size subsets. Default = [10, 90].')
    parser.add_argument('-amp', '--amplifier', default=5, type=float,
        help='Amplifier value to plot. Default = 5.')
    parser.add_argument('-do', '--ADO', nargs='+', default=[0, 0.2, 0.4],
        type=float,
        help='Simulated ADO value to plot. Default = [0, 0.2, 0.4].')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        default=['cellcoal', 'cellphy'],
        help='Method to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_PRC(args)