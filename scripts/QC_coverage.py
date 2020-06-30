#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
from itertools import cycle
from matplotlib import pyplot as plt


COLORS = cycle([
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', #light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
])


def parse_args():
    parser = argparse.ArgumentParser(
        prog='QC_coverage', usage='python3 QC_coverage.py <DATA> [options]',
        description='*** Generate Lorenz curve and Gini coefficient. ***'
    )
    parser.add_argument(
        'input', nargs='+', type=str,
        help='Absolute or relative path(s) to input data in BED format'
    )
    parser.add_argument(
        '-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.'
    )
    
    args = parser.parse_args()
    return args


def calc_gini(x):
    return -1
    n = x.shape[0]
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n


def main(args):
    sample_names = [os.path.basename(i).split('.')[0] for i in args.input]
    cols = ['Total Depth', 'Cov. Depth', 'Total Breadth', '5x Total Breadth',
        '5x Cov. Breadth', '10x Total Breadth', '10x Cov. Breadth',
        '20x Total Breadth', '20x Cov. Breadth', 'Gini Coefficient']
    df_sum = pd.DataFrame(index=sample_names + ['mean', 'stdev'], columns=cols)
    df_sum.index.name = 'Sample'

    fig1, ax1 = plt.subplots(figsize=[6,6])
    fig2, ax2 = plt.subplots(figsize=[6,6])

    for i, sample_file in enumerate(args.input):
        sample_name = sample_names[i]

        # 0: genome/chromosom, 1: coverage depth, 2: # bases with covg 1,
        # 3: size of total genome/chromosome, 4: fraction of bases with covg 1
        df = pd.read_csv(sample_file, delimiter='\t', header=None)
        # Add column with: coverage [reads] per bp [depth]
        df[df.shape[1] + 1] = df[1] * df[2]
        
        avg_depth = np.average(df[1].values, weights=df[2].values)
        cov_depth = np.average(df.iloc[1:,1].values, weights=df.iloc[1:,2].values)
        cov_total = 1 - df.iloc[0, 4]

        # normalized over covered basepairs
        x = np.cumsum((df.iloc[1:,4] / cov_total).values)
        x = np.insert(x, 0, 0) 
        # normalized over full genome/exome
        # x = np.cumsum(df[4].values)
        
        bp_total = df[6].sum()
        y = np.cumsum(df[6].values) / bp_total

        cov_breadth = []
        for cov_x in [5, 10, 20]:
            breadth = df.iloc[cov_x:, 4].sum()
            cov_breadth.extend([breadth, breadth / cov_total])

        df_sum.iloc[i] = [avg_depth, cov_depth, cov_total] + cov_breadth \
            + [calc_gini(x)]

        sample_col = next(COLORS)
        # Scatter plot of Lorenz curve
        ax1.plot(x, y, color=sample_col, label=sample_name)
        # Read depth distribution
        ax2.plot(df.iloc[1:, 1].values, (1 - x)[:-1], color=sample_col,
            label=sample_name)
        
    df_sum.loc['mean'] = df_sum.mean()
    df_sum.loc['stdev'] = df_sum.std()

    if args.output != '':
        out_dir = args.output
    else:
        out_dir = os.path.dirname(sample_file)

    if df.iloc[0, 3] > 3e9:
        seq = 'Exome'
    else:
        seq = 'Genome'
    plot_lorenz(fig1, ax1, out_dir, x_label=seq)
    plot_read_depth_dist(fig2, ax2, out_dir, x_max=(x < 0.99).sum())
    plt.close()

    df_sum.to_csv(os.path.join(out_dir, 'QC_sequencing.tsv'), sep='\t')


def plot_lorenz(fig, ax, out_dir, x_label='Genome'):
    ax.plot([0, 1], [0, 1], '--', color='k', label='Perfect Uniformity')
    ax.set_xlabel(f'Cumulative Fraction of Covered {x_label}')
    ax.set_ylabel('Cumulative Fraction of Total Basepairs')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.legend()
    ax.grid(True)
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9)

    fig.savefig(os.path.join(out_dir, 'Lorenz.pdf'), dpi=300)
    

def plot_read_depth_dist(fig, ax, out_dir, x_max=None):
    ax.set_xlabel('Coverage Depth')
    ax.set_ylabel('Fraction of Coverage $\geq$ depth')
    ax.set_xlim([0, x_max])
    ax.set_ylim([0, 1])
    ax.legend()
    ax.grid(True)
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9)

    fig.savefig(os.path.join(out_dir, 'ReadDepthDist.pdf'), dpi=300)


if __name__ == '__main__':
    args = parse_args()
    main(args)