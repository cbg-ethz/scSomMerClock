#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from scipy.stats import chi2

from defaults import *


def generate_pval_plot_clock(args):
    df = pd.DataFrame(columns=['ADO', 'wMax', 'method', 'P-value', 'Lambda'])

    for res_file in sorted(os.listdir(args.input)):
        if not res_file.startswith('res_clock0') or 'bulk' in res_file:
            continue

        if re.match('res_clock\d+x_', res_file):
            continue

        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO not in args.ADO:
            continue
        print(f'Including file: {res_file}')

        df_summary = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_summary.drop([-1], inplace=True)
        df_summary = df_summary.T.reset_index().drop_duplicates() \
            .set_index('index').T

        cols = [i for i in df_summary.columns if i.startswith('p-value') \
            and not 'poissonTree_0' in i]

        for col in cols:
            if col[8:].startswith('poissonTree'):
                try:
                    wMax = int(re.search('wMax(\d+)', col).group(1))
                except AttributeError:
                    # No errors simulated
                    wMax = args.wMax
                else:
                    if ADO == 0 and len(cols) < 6:
                        wMax = args.wMax
                    else:
                        if wMax not in args.wMax:
                            continue
                        wMax = [wMax]
                method = col.split('.')[-1]

                if method not in COLORS:
                    method = col.split('_')[-1]

                if method not in args.method:
                    continue
            elif col[8:].startswith('poissonDisp'):
                if 'poissonDisp' not in args.method:
                    continue
                wMax = [-1]
                method = 'poissonDisp'
            elif col[8:].startswith('paup'):
                if 'PAUP*' not in args.method:
                    continue
                wMax = [-1]
                method = col.split('.')[-1]
                if method.startswith('p-value'):
                    method = col.split('_')[-1]
                if method == 'cellcoal':
                    method = 'PAUP*'
                else:
                    continue
            else:
                raise TypeError(f'Unknown column type: {col}')

            df_new = df_summary.loc[:,[col, col.replace('p-value_', '-2logLR_')]]
            df_new.columns = ['P-value', 'Lambda']

            df_new['ADO'] = ADO
            df_new['method'] = method

            for wMax_i in wMax:
                df_new['wMax'] = wMax_i
                df = df.append(df_new, ignore_index=True)

    for col in ['ADO', 'wMax', 'P-value', 'Lambda']:
        df[col] = df[col].astype(float)

    wMax_vals = df['wMax'].unique()
    wMax_vals.sort()
    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()

    if args.statistic:
        plt_fct = plot_lambda_dist
    else:
        plt_fct = plot_pVal_dist

    row_no = wMax_vals.size
    col_no = ADO_vals.size
    fig, axes = get_subplots(row_no, col_no)

    for i, ADO_val in enumerate(ADO_vals):
        df_plot = df[df['ADO'] == ADO_val]
        plt_fct(df_plot, wMax_vals, axes[:,i], (i, col_no))

        if col_no > 1:
            if ADO_val > 0:
                col_title = f'FN rate: {(ADO_val / 2) * 100: >4.1f}%'
            else:
                col_title = 'No Errors'
            add_col_header(axes[0, i], col_title)

    if wMax_vals.size > 2:
        for j, wMax in enumerate(wMax_vals):
            if wMax == -1:
                continue
            add_row_header(axes[j, -1], r'$w_{{max}}=$' + f'\n{wMax:.0f}')

    plot_fig(fig, args.output)


def plot_pVal_dist(df, wMaxs, axes, col_no):
    for i, wMax in enumerate(wMaxs):
        ax = axes[i]
        data = df[df['wMax'] == wMax]
        if wMax == -1:
            if not 'PAUP*' in data['method'].unique():
                sig_data = np.repeat(
                    [[data['ADO'].unique()[0], -1, 'PAUP*', 0.001, 99]],
                    data.shape[0],
                    axis=0
                )
                dtypes = data.dtypes.to_dict()
                data = data.append(pd.DataFrame(sig_data, columns=data.columns))
                data = data.astype(dtypes)

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

        # ax.annotate(f'n = {data["method"].value_counts().min():.0f}',
        #     xy=(0.9, 0.825), xytext=(0, 5),  xycoords='axes fraction',
        #     textcoords='offset points', ha='right', va='top', bbox=bbox_props)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Not last row
        if i < wMaxs.size - 1 and wMax != -1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        # First col
        if col_no[0] == 0:
            ax.set_ylabel('Probability')
        else:
            ax.set_ylabel(None)
            ax.set_yticklabels([])

        # Middle col
        if col_no[0] == np.floor(col_no[1] / 2) \
                and (i == wMaxs.size - 1 or wMax == -1):
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')

        # Last col
        # if col_no[0] == col_no[1] - 1:
        #     ax2 = ax.twinx()
        #     ax2.set_ylabel('Clock')
        #     ax2.set_yticks([])
        #     for tick in  ax.yaxis.majorTicks:
        #         tick.tick1line.set_markersize(0)


def plot_lambda_dist(df, wMax, axes, col_no):
    dof = 29 # TODO <NB> remove hardcoding
    bins = 100

    x = np.linspace(chi2.ppf(0.001, dof), chi2.ppf(0.999, dof), bins)
    y = chi2.pdf(x, dof)


    for i, j in enumerate(wMax):
        ax = axes[i]
        data = df[df['wMax'] == j]

        dp = sns.histplot(data, x='Lambda', hue='method', stat='density',
            bins=bins, common_norm=False, common_bins=False, legend=False,
            hue_order=HUE_ORDER, palette=COLORS, ax=ax,
        )
        ax.plot(x, y, ls='--', color='black', label=f'Chi2({dof}) pdf')

        ax.set_xlim((0, 100))
        ax.set_ylim((0, 0.2))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < wMax.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        # First col
        if col_no[0] == 0:
            if i != np.floor(wMax.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel(f'Density')
        else:
            ax.set_ylabel(None)
            ax.set_yticklabels([])

        # Middle col
        if col_no[0] == np.floor(col_no[1] / 2) and i == wMax.size - 1:
            ax.set_xlabel(r'$\lambda$')
        else:
            ax.set_xlabel('')

        # Last col
        # if col_no[0] == col_no[1] - 1:
        #     ax2 = ax.twinx()
        #     ax2.set_ylabel('Clock')
        #     ax2.set_yticks([])
        #     for tick in  ax.yaxis.majorTicks:
        #         tick.tick1line.set_markersize(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-spd', '--skip_poisson', action='store_true',
        help='Skip Poisson distribution test.')
    parser.add_argument('-w', '--wMax', nargs='+', type=float, default=[1000],
        help='wMax values to plot. Default = 1000.')
    parser.add_argument('-a', '--ADO', nargs='+', type=float,
        default=[0, 0.05, 0.1, 0.2, 0.4, 0.6],
        help='ADO values to plot. Default = [0, 0.05, 0.1, 0.2, 0.4, 0.6].')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy', 'scite', 'poissonDisp', 'PAUP*'],
        default=['cellcoal', 'cellphy', 'poissonDisp', 'PAUP*'],
        help='Method to plot. Default = all.')
    parser.add_argument('-s', '--statistic', action='store_true',
        help='Plot test statistic as separate figure.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_clock(args)