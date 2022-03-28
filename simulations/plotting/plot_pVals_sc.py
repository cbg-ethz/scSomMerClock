#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import chi2

from defaults import *


def generate_pval_plot(args):
    df = pd.DataFrame(columns=['ADO', 'Amplifier', 'Tree', 'P-value', 'Lambda'])

    for res_file in os.listdir(args.input):
        if not res_file.startswith('res_clock') or 'bulk' in res_file:
            continue
        ADO = float(re.search('WGA(0[\.\d]*)', res_file).group(1))

        if ADO not in args.ADO:
            continue

        if res_file.startswith(('res_clock00', 'res_clock0_')):
            ampl = 1
        else:
            ampl = float(re.search('res_clock(\d[\.\d]*)x', res_file).group(1))
        if args.amplifier and ampl not in args.amplifier:
            continue

        df_summary = pd.read_csv(
            os.path.join(args.input, res_file), sep='\t', index_col='run')
        df_summary.drop([-1], inplace=True)
        df_summary = df_summary.T.reset_index().drop_duplicates() \
            .set_index('index').T

        cols = [i for i in df_summary.columns if i.startswith('p-value')]
        for col in cols:
            if col[8:].startswith('poissonTree'):
                if ADO == 0 and ampl == 1:
                    if col.split('.')[0][-1] == '0':
                        continue
                    else:
                        pass
                else:
                    try:
                        wMax = int(re.search('wMax(\d+)', col).group(1))
                    except AttributeError:
                        continue
                    if wMax != args.wMax:
                        continue
                tree = col.split('.')[-1]
                if tree not in colors:
                    tree = col.split('_')[-1]
                if tree not in args.method:
                    continue
            elif col[8:].startswith('poissonDisp'):
                if 'poissonDisp' not in args.method or ampl > 1:
                    continue
                tree = '-'
            elif col[8:].startswith('paup'):
                if 'PAUP*' not in args.method:
                    continue
                tree = col.split('.')[-1]
                if tree == 'cellcoal':
                    tree = 'PAUP*'
                else:
                    continue
            else:
                raise TypeError(f'Unknown column type: {col}')

            if ampl > 1:
                df_summary = df_summary[
                    (df_summary['aff. cells'] > 2) \
                        & (df_summary['aff. cells'] < 28)]
            df_new = df_summary.loc[:,[col, col.replace('p-value_', '-2logLR_')]]
            df_new.columns = ['P-value', 'Lambda']
            df_new['ADO'] = ADO
            df_new['Amplifier'] = ampl
            df_new['Tree'] = tree

            df = df.append(df_new, ignore_index=True)

    ampl_vals = df['Amplifier'].unique()
    ampl_vals.sort()
    ADO_vals = df['ADO'].unique()
    ADO_vals.sort()

    if args.legend:
        generate_legend_plot(df['Tree'].unique(), args.output)

    fig, axes = plt.subplots(nrows=ampl_vals.size, ncols=ADO_vals.size,
        figsize=(3 * ADO_vals.size, ampl_vals.size + 1))
    fig2, axes2 = plt.subplots(nrows=ampl_vals.size, ncols=ADO_vals.size,
        figsize=(3 * ADO_vals.size, ampl_vals.size + 1))
    axes = np.reshape(axes, (ampl_vals.size, ADO_vals.size))
    axes2 = np.reshape(axes2, (ampl_vals.size, ADO_vals.size))

    for i, ADO_val in enumerate(ADO_vals):
        df_plot = df[df['ADO'] == ADO_val]
        plot_pVal_dist(df_plot, ampl_vals, axes[:,i], (i, ADO_vals.size))
        plot_lambda_dist(df_plot, ampl_vals, axes2[:,i], (i, ADO_vals.size))

        if ADO_val == 0:
            col_title = 'No Errors'
        else:
            col_title = f'FN rate: {ADO_val / 2}'

        axes[0,i].annotate(col_title, xy=(0.5, 1.15), xytext=(0, 5),
            xycoords='axes fraction', textcoords='offset points', size='large',
            ha='center', va='baseline')
        axes2[0,i].annotate(col_title, xy=(0.5, 1.15), xytext=(0, 5),
            xycoords='axes fraction', textcoords='offset points', size='large',
            ha='center', va='baseline')

    fig.tight_layout()
    fig2.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=DPI)
        fig2.savefig(os.path.splitext(args.output)[0] + '_lambda.png', dpi=DPI)
    else:
        plt.show()
    plt.close()



def plot_pVal_dist(df, ampl_vals, axes, col):
    for i, ampl in enumerate(ampl_vals):
        ax = axes[i]
        data = df[df['Amplifier'] == ampl]

        hue_order = sorted(data['Tree'].unique())
        if len(hue_order) > 1:
            hue_order.append('gap')

        dp = sns.histplot(data, x='P-value', hue='Tree',
            element='bars', stat='probability', kde=False,
            common_norm=False, fill=True,
            binwidth=0.05, binrange=(0, 1), multiple='layer',
            palette=colors,
            hue_order=hue_order,
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
            add_rugs(rug_data, offset=k, ax=ax, color=colors[method])
            k += 1

        ax.spines['top'].set_visible(False)

        n = (data['Tree'] == 'cellcoal').sum()
        ax.annotate(f'n = {n}', xy=(0.9, 0.75), xytext=(0, 5),
            xycoords='axes fraction', textcoords='offset points',
            ha='right', va='top', bbox=bbox_props)

        if i < ampl_vals.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

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
                ax2.set_ylabel(f'\nAmplifier:\n{ampl:.0f}x')
            ax2.set_yticks([])


def plot_lambda_dist(df, ampl_vals, axes, col):
    dof = 29 # TODO <NB> remove hardcoding
    bins = 100

    x = np.linspace(chi2.ppf(0.001, dof), chi2.ppf(0.999, dof), bins)
    y = chi2.pdf(x, dof)


    for i, ampl in enumerate(ampl_vals):
        ax = axes[i]
        data = df[df['Amplifier'] == ampl]

        dp = sns.histplot(data, x='Lambda', hue='Tree', stat='density', bins=bins,
            common_norm=False, common_bins=False, legend=False,
            hue_order=HUE_ORDER, palette=colors, ax=ax,
        )
        ax.plot(x, y, ls='--', color='black', label=f'Chi2({dof}) pdf')

        ax.set_xlim((0, 100))
        ax.set_ylim((0, 0.2))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if i < ampl_vals.size - 1:
            ax.set_xticklabels([])
            ax.set_xlabel(None)

        if col[0] == np.floor((col[1] - 1) / 2) and i == len(ampl_vals) - 1:
            ax.set_xlabel(r'$\lambda$')
        else:
            ax.set_xlabel('')

        if col[0] == 0 and i == np.floor((len(ampl_vals) - 1) / 2):
            ax.set_ylabel('Density')
        else:
            ax.set_ylabel('')

        if col[0] == col[1] - 1:
            ax2 = ax.twinx()
            if ampl == 1:
                ax2.set_ylabel(f'\nClock')
            else:
                ax2.set_ylabel(f'\nAmplifier:\n{ampl:.0f}x')
            ax2.set_yticks([])


def generate_legend_plot(methods, output):
    fig, ax = plt.subplots(figsize=(4, 3))

    handles = []
    labels = []
    for method in HUE_ORDER:
        if method not in methods:
            continue
        labels.append(methods_names[method])
        #handles.append(Line2D([0], [0], color=colors[method], lw=2))
        handles.append(plt.Rectangle((0,0),1,1, color=colors[method]))

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
    parser.add_argument('-w', '--wMax', type=float, default=400,
        help='wMax value to plot. Default = all.')
    parser.add_argument('-a', '--ADO', nargs='+', type=float,
        default=[0, 0.2, 0.4, 0.6],
        help='ADO values to plot (columns). Default = all.')
    parser.add_argument('-amp', '--amplifier', nargs='+', type=float,
        default=[1, 2, 5],
        help='Amplifier values to plot. Clock = 1. Default = [1, 5, 10].')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy', 'scite', 'poissonDisp', 'PAUP*'],
        default=['cellcoal', 'scite', 'cellphy', 'poissonDisp', 'PAUP*'],
        help='Method to plot. Default = all.')
    parser.add_argument('-l', '--legend', action='store_true',
        help='Plot legend as separate figure.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_pval_plot(args)