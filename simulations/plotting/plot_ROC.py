#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

from defaults import *


def generate_PRC(args):
    cols =  ['ADO', 'aff. cells', 'method', 'wMax', 'predicted']
    df = pd.DataFrame(columns=cols)

    for res_file in sorted(os.listdir(args.input)):
        if not res_file.startswith('res_clock') or 'bulk' in res_file:
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO not in args.ADO:
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
                if method not in args.method:
                    continue
                if ADO == 0 and clock:
                    wMax = int(rel_col.split('_')[-1].split('.')[0])
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

    min_cells = np.array(args.min_cell)
    ADO_vals = df['ADO'].unique()

    single_plot = min_cells.size == 1 and ADO_vals.size == 1

    if single_plot:
        fig1, axes1 = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
        fig2, axes2 = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig1, axes1 = plt.subplots(nrows=min_cells.size, ncols=ADO_vals.size,
            figsize=(3 * ADO_vals.size, min_cells.size * 3))
        fig2, axes2 = plt.subplots(nrows=min_cells.size, ncols=ADO_vals.size,
            figsize=(3 * ADO_vals.size, min_cells.size * 3))
    axes1 = np.reshape(axes1, (min_cells.size, ADO_vals.size))
    axes2 = np.reshape(axes2, (min_cells.size, ADO_vals.size))

    for i, ADO_val in enumerate(ADO_vals):
        df_plot = df[df['ADO'] == ADO_val]
        plot_curves(df_plot, axes1, axes2, min_cells, (i, ADO_vals.size), args.amplifier)

        if not single_plot:
            if ADO_val == 0:
                col_title = 'No Errors'
            else:
                col_title = f'FN rate: {ADO_val / 2}'

            axes1[0,i].annotate(col_title, xy=(0.5, 1.15), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points', size='large',
                ha='center', va='baseline')
            axes2[0,i].annotate(col_title, xy=(0.5, 1.15), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points', size='large',
                ha='center', va='baseline')

    fig1.tight_layout()
    fig2.tight_layout()

    if args.output:
        fig1.savefig(f'PRC_{args.output}.png', dpi=300)
        fig2.savefig(f'ROC_{args.output}.png', dpi=300)
    else:
        plt.show()
    plt.close()


def plot_curves(df_in, axes1, axes2, min_cells, col_id, ampl):
    wMax_map = {j: i for i, j in enumerate(sorted(df_in['wMax'].unique()))}
    cmap = cm.get_cmap('viridis_r', len(wMax_map))

    for j, min_cell in enumerate(min_cells):
        df_cell = df_in[(df_in['aff. cells'] >= min_cell) \
        & ((df_in['aff. cells'] <= 30 - min_cell) | (df_in['aff. cells'] == np.inf))]

        if df_cell.size == 0:
            print(f'!WARNING - No run with {min_cell} cells affected!')
            continue

        ax1 = axes1[j, col_id[0]]
        ax2 = axes2[j, col_id[0]]

        for method, df in df_cell.groupby('method'):
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
                    fprs.append(fpr)
            else:
                vals = df['predicted'].value_counts()
                rec = [vals.TP / (vals.TP + vals.FN)]
                prec = [vals.TP / (vals.TP + vals.FP)]

            sns.lineplot(x=rec, y=prec, ax=ax1, color=colors[method],
                markersize=0, alpha=0.75) #PRC
            sns.lineplot(x=fprs, y=rec, ax=ax2, color=colors[method],
               markersize=0, alpha=0.75) # ROC

            for k, wMax in enumerate(w_max, 2):
                ax1.plot(rec[k], prec[k], marker='o', markersize=6,
                    mec=colors[method], mfc=cmap(k - 2))
                ax2.plot(fprs[k], rec[k], marker='o', markersize=6,
                    mec=colors[method], mfc=cmap(k - 2))
            ax2.plot([0, 1], [0, 1], color='grey', ls='--')

        ax1.set_xlim((-0.05, 1.05))
        ax1.set_ylim((0.4, 1.05))

        ax2.set_xlim((-0.05, 1.05))
        ax2.set_ylim((-0.05, 1.05))

        # First col, middle row
        if col_id[0] == 0 and j == np.floor(min_cells.size / 2):
            ax1.set_ylabel('Precision')
            ax2.set_ylabel('True positive rate')
        else:
            ax1.set_ylabel('')
            ax2.set_ylabel('')

        # Middle col, last row
        if col_id[0] == np.floor((col_id[1] - 1) / 2) and j == min_cells.size - 1:
            ax1.set_xlabel('Recall')
            ax2.set_xlabel('False positive rate')
        else:
            ax1.set_xlabel('')
            ax2.set_xlabel('')

        # Last column
        if col_id[0] == col_id[1] - 1:
            ax12 = ax1.twinx()
            ax12.set_yticks([])
            ax22 = ax2.twinx()
            ax22.set_yticks([])
            if len(min_cells) == 1:
                ax12.set_ylabel(f'\nAmplifier:\n{ampl:.0f}', fontsize=8)
                ax22.set_ylabel(f'\nAmplifier:\n{ampl:.0f}', fontsize=8)
            else:
                ax12.set_ylabel('\n' + r'$\geq$' + f'{min_cell} cells', fontsize=8)
                ax22.set_ylabel('\n' + r'$\geq$' + f'{min_cell} cells', fontsize=8)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Directory with summary files.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-c', '--min_cell', nargs='+', default = [1, 2, 3, 4, 5],
        type=float, help='Min. #cells affected. Default = [1, 2, 3, 4, 5].')
    parser.add_argument('-a', '--amplifier', default=2, type=float,
        help='Amplifier value to plot. Default = 2.')
    parser.add_argument('-do', '--ADO', nargs='+', default=[0, 0.2, 0.4, 0.6],
        type=float,
        help='Simulated ADO value to plot. Default = [0, 0.2, 0.4, 0.6].')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        default=['cellcoal', 'scite', 'cellphy'],
        help='Method to plot. Default = all.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_PRC(args)