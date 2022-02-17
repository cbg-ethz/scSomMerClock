#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from scipy.stats import chi2

from defaults import *


def generate_pval_plot_clock(args):
    df = pd.DataFrame(columns=['ADO', 'wMax', 'Tree', 'P-value', 'Lambda'])

    for res_file in os.listdir(args.input):
        if not res_file.startswith('res_clock0') or 'bulk' in res_file:
            continue
        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))
        if ADO not in args.ADO:
            continue

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
            elif col[8:].startswith('poissonDisp'):
                if args.skip_poisson:
                    continue
                wMax = [-1]
                tree = '-'
            elif col[8:].startswith('paup'):
                if args.skip_poisson:
                    continue
                wMax = [-1]
                tree = 'PAUP*'
            else:
                raise TypeError(f'Unknown column type: {col}')

            df_new = df_summary.loc[:,[col, col.replace('p-value_', '-2logLR_')]]
            df_new.columns = ['P-value', 'Lambda']

            df_new['ADO'] = ADO
            df_new['Tree'] = tree
            for wMax_i in wMax:
                df_new['wMax'] = wMax_i
                df = df.append(df_new, ignore_index=True)

    # If single ADO and wMax value, plot all available data in one plot
    single_plot = len(args.wMax) == 1 and len(args.ADO) == 1
    if single_plot:
        df.loc[df['wMax'] == -1, 'wMax'] = args.wMax[0]

    wMax_vals = df['wMax'].unique()
    wMax_vals.sort()
    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()

    if args.legend:
        generate_legend_plot(df['Tree'].unique(), args.output)

    if single_plot:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
        fig2, axes2 = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    else:
        fig, axes = plt.subplots(nrows=wMax_vals.size, ncols=ADO_vals.size,
            figsize=(3 * ADO_vals.size, wMax_vals.size + 1))
        fig2, axes2 = plt.subplots(nrows=wMax_vals.size, ncols=ADO_vals.size,
            figsize=(3 * ADO_vals.size, wMax_vals.size + 1))
    axes = np.reshape(axes, (wMax_vals.size, ADO_vals.size))
    axes2 = np.reshape(axes2, (wMax_vals.size, ADO_vals.size))

    for i, ADO_val in enumerate(ADO_vals):
        if single_plot:
            col_type = 'only'
        elif i == 0:
            col_type = 'first'
        elif i == ADO_vals.size - 1:
            col_type = 'last'
        elif i == np.floor(ADO_vals.size / 2):
            col_type = 'middle'
        else:
            col_type = 'intermediate'

        df_plot = df[df['ADO'] == ADO_val]
        plot_pVal_dist(df_plot, wMax_vals, axes[:,i], col_type)
        plot_lambda_dist(df_plot, wMax_vals, axes2[:,i], col_type)

        if not single_plot:
            axes[0,i].annotate(f'ADO rate: {ADO_val}', xy=(0.5, 1.1),
                xytext=(0, 5), xycoords='axes fraction',
                textcoords='offset points', size='large', ha='center',
                va='baseline')
            axes2[0,i].annotate(f'ADO rate: {ADO_val}', xy=(0.5, 1.1),
                xytext=(0, 5), xycoords='axes fraction',
                textcoords='offset points', size='large', ha='center',
                va='baseline')

    if single_plot:
        fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.75)
        fig2.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    else:
        fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
            hspace=0.5, wspace=0.5)
        fig2.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
            hspace=0.5, wspace=0.5)

    if args.output:
        fig.savefig(args.output, dpi=DPI)
        fig2.savefig(os.path.splitext(args.output)[0] + '_lambda.png', dpi=DPI)
    else:
        plt.show()
    plt.close()



def plot_pVal_dist(df, wMax, axes, column_type):
    for i, j in enumerate(wMax):
        ax = axes[i]
        data = df[df['wMax'] == j]
        dp = sns.histplot(data, x='P-value', hue='Tree',
            element='poly', stat='probability', kde=False,
            common_norm=False, fill=False,
            binwidth=0.01, binrange=(0, 1), multiple='dodge',
            kde_kws={'cut': 0, 'clip': (0, 1)},
            line_kws={'lw': 3, 'alpha': 0.66},
            palette=colors, alpha=0.5,
            hue_order=HUE_ORDER,
            legend=False, ax=ax,
        )

        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))

        # Add rugplots
        k = 0
        for method in HUE_ORDER:
            rug_data = data[data['Tree'] == method]['P-value'].values
            if rug_data.size == 0:
                continue

            segs = np.stack((np.c_[rug_data, rug_data],
                np.c_[np.zeros_like(rug_data) + 1 + RUG_HEIGHT*2 * k,
                        np.zeros_like(rug_data) + 1 + RUG_HEIGHT*2 * (k + 1)]),
                    axis=-1)
            lc = LineCollection(segs, transform=ax.get_xaxis_transform(),
                clip_on=False, color=colors[method], linewidth=0.05, alpha=0.75)
            ax.add_collection(lc)
            k += 1

        l = 0
        for tree, col in colors.items():
            df_sub = df[(df['Tree'] == tree) & (df['wMax'] == j)]
            if df_sub.size > 0:
                y_pval = (df_sub['P-value'] <= 0.05).mean()
                if y_pval > 0.33:
                    va = 'top'
                    y_dist = -0.02
                else:
                    va = 'bottom'
                    y_dist = 0.02

                ax.axhline(y_pval, ls='--', color=col, lw=1)
                # ax.axvline(0.05, ls='--', color='grey', lw=1)
                ax.text(0.05 + l * 0.15, y_pval + y_dist, f'{y_pval:.2f}',
                    color=col, ha='left', va=va, rotation=45)
                l += 1

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < wMax.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if column_type == 'first':
            if i != np.floor(wMax.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel('Probability')
        else:
            ax.set_ylabel('')

        if column_type == 'middle' and i == wMax.size - 1:
            ax.set_xlabel('P-value')
        else:
            ax.set_xlabel('')

        if column_type == 'last':
            ax2 = ax.twinx()
            if j >= 0:
                ax2.set_ylabel('\n' + r'$w_{max}=$ '+ f'\n{j:.0f}', fontsize=12)
            else:
                ax2.set_ylabel('\nPoisson\nDispersion', fontsize=12)
            ax2.set_yticks([])

        if column_type == 'only':
            ax.set_ylabel('Probability')
            ax.set_xlabel('P-value')


def plot_lambda_dist(df, wMax, axes, column_type):
    dof = 29 # TODO <NB> remove hardcoding
    bins = 100

    x = np.linspace(chi2.ppf(0.001, dof), chi2.ppf(0.999, dof), bins)
    y = chi2.pdf(x, dof)


    for i, j in enumerate(wMax):
        ax = axes[i]
        data = df[df['wMax'] == j]

        dp = sns.histplot(data, x='Lambda', hue='Tree', stat='density', bins=bins,
            common_norm=False, common_bins=False, legend=False,
            hue_order=HUE_ORDER, palette=colors, ax=ax,
        )
        ax.plot(x, y, ls='--', color='black', label=f'Chi2({dof}) pdf')

        ax.set_xlim((0, 100))
        ax.set_ylim((0, 0.2))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < wMax.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if column_type == 'first':
            if i != np.floor(wMax.size / 2):
                ax.set_ylabel('')
            else:
                ax.set_ylabel(f'Density')
        else:
            ax.set_ylabel('')

        if column_type == 'middle' and i == wMax.size - 1:
            ax.set_xlabel(r'$\lambda$')
        else:
            ax.set_xlabel('')

        if column_type == 'last':
            ax2 = ax.twinx()
            if j >= 0:
                ax2.set_ylabel('\n' + r'$w_{max}=$ '+ f'\n{j:.0f}', fontsize=12)
            else:
                ax2.set_ylabel('\nPoisson\nDispersion', fontsize=12)
            ax2.set_yticks([])

        if column_type == 'only':
            ax.set_xlabel(r'$\lambda$')
            ax.set_ylabel(f'Density')


def generate_legend_plot(methods, output):
    fig, ax = plt.subplots(figsize=(4, 3))

    handles = []
    labels = []
    for method in HUE_ORDER:
        if method not in methods:
            continue
        labels.append(methods_names[method])
        handles.append(Line2D([0], [0], color=colors[method], lw=2))
    handles.extend([Line2D([0], [0], color='white'),
        Line2D([0], [0], color='black', lw=2, ls='--')])
    labels.extend(['', r'% of p-values $\leq$ 0.05'])

    ax.grid(False)
    ax.axis('off')
    ax.legend(handles, labels, frameon=True, title=r'$\bf{Method}$')

    fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    if output:
        fig.savefig(os.path.splitext(output)[0] + '_legend.png', dpi=DPI)
    else:
        plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-spd', '--skip_poisson', action='store_true',
        help='Skip Poisson distribution test.')
    parser.add_argument('-w', '--wMax', nargs='+', type=float,
        default=[1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
        help='wMax values to plot. Default = all.')
    parser.add_argument('-a', '--ADO', nargs='+', type=float,
        default=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
        help='ADO values to plot. Default = all.')
    parser.add_argument('-l', '--legend', action='store_true',
        help='Plot legend as separate figure.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot_clock(args)