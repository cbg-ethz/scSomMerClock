#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd
from scipy.stats import pearsonr

from defaults import *


def generate_scatter_plot(args):
    df = pd.read_csv(args.input, sep='\t', na_values=['?', 'NaN'])
    df.dropna(axis=1, how='all', inplace=True)
    df.drop(['FP', 'Signatures'], axis=1, inplace=True)
    df.rename({'cells': '# cells', 'Poisson Tree Test': 'PTT [p-val]',
        'Clock-like Signatures [%]': 'SBS 1 & SBS5 [%]'}, inplace=True,axis=1)
    df['subset'] = df['subset'].str.replace('Neurons', 'Normal') \
        .replace('all', 'Normal')

    parameters = df.columns.values[2:]
    df[parameters] = df[parameters].astype(float)

    row_no = parameters.size
    col_no = parameters.size
    fig, axes = get_subplots(row_no, col_no)

    for i, par1 in enumerate(parameters):

        for j, par2 in enumerate(parameters):
            df_plot = df[['dataset', 'subset', par2, par1]]
            plot_scatter(df_plot, axes[i, j], (i, row_no), (j, col_no),
                i == 0 and j == col_no - 1)

            if i == 0:
                add_col_header(axes[0, j], par2)
                add_row_header(axes[j, -1], par2)

    MARGINS = {
        'left': 0.15,
        'right': 0.7,
        'top': 0.8,
        'bottom': 0.3,
        'wspace': 0.5,
        'hspace': 0.75,
    }
    plt.subplots_adjust(**MARGINS)
    if args.output:
        if not args.output.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            args.output += '.png'
        fig.savefig(args.output, dpi=DPI)
    else:
        plt.show()


def plot_scatter(df, ax, row_ids, col_ids, legend=False):
    x_col = df.columns[2]
    y_col = df.columns[3]

    if x_col != y_col:
        scp = sns.scatterplot(data=df, x=x_col, y=y_col, style='subset',
            hue='dataset', ax=ax, legend=legend)

        df = df.dropna()
        corr = pearsonr(df[x_col].values, df[y_col].values)
        ax.annotate(f'PCC={corr[0]:.2f}', xy=(0.9, 0.825), xytext=(0, 5),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='top')

    if legend:
        l_handels, l_labels = ax.get_legend_handles_labels()
        ax.legend(l_handels, l_labels, frameon=True,
            bbox_to_anchor=(1.4, 1), loc=2, borderaxespad=0.)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if col_ids[0] == 0:
        ax.set_ylabel(y_col)
    else:
        ax.set_ylabel('')

    if row_ids[0] == row_ids[1] - 1:
        ax.set_xlabel(x_col)
    else:
        ax.set_xlabel('')



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    generate_scatter_plot(args)