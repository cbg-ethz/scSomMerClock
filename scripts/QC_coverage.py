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
    cols = ['Avg. Depth', 'Cov. Depth', 'Avg. Breadth',
        '5x Cov. Breadth', '10x Cov. Breadth', '20x Cov. Breadth',
        'Gini Coefficient']
    df_sum = pd.DataFrame(index=sample_names + ['average', 'stdev'], columns=cols)
    df_sum.index.name = 'Sample'

    fig1, ax1 = plt.subplots(figsize=[6,6])
    fig2, ax2 = plt.subplots(figsize=[6,6])

    for i, sample_file in enumerate(args.input):
        sample_name = sample_names[i]

        data = np.genfromtxt(sample_file, delimiter='\t')
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
        # normalized over full genome/exome
        # x = np.cumsum(df[4].values)
        x = np.insert(x, 0, 0) 

        bp_total = df[6].sum()
        y = np.cumsum(df[6].values) / bp_total
        y = np.insert(y, 0, 0)
        import pdb; pdb.set_trace()

        cov_breadth = []
        for cov_x in [5, 10, 20]:
            cov_breadth.append(
                np.where(data[:, 1] >= cov_x, data[:,4], 0).sum() / cov_total
            )
            # import pdb; pdb.set_trace()
            test = np.average(data[cov_x:, 1], weights=data[cov_x:,2])
            # print(cov_x, test)

        df_sum.iloc[i] = [avg_depth, cov_depth, cov_total] + cov_breadth \
            + [calc_gini(x)]

        sample_col = next(COLORS)
        # Scatter plot of Lorenz curve
        ax1.plot(x, y, color=sample_col, label=sample_name)
        # Read depth distribution
        ax2.plot(data[1:, 1], (1 - x)[:-1], color=sample_col, label=sample_name)
        # Summary with: sample name, Avg depth, 1/5/10/20x breadth, gini coeff
        
    import pdb; pdb.set_trace()

    if args.output != '':
        out_dir = args.output
    else:
        out_dir = os.path.dirname(sample_file)
    out_name = sample_name.split('_')[0]

    plot_lorenz(fig1, ax1, out_dir, out_name)
    plot_read_depth_dist(fig2, ax2, out_dir, out_name)
    stdout_summary(summary, out_dir, out_name)
    plt.close()

def stdout_summary(summary, out_dir, out_name, sep='\t'):
    out_file = os.path.join(out_dir, f'{out_name}_QC_sequencing.tsv')
    with open(out_file, 'w') as f:
        f.write(sep.join(['Sample', 'Avg. Depth', 'Cov. Depth',
            'Avg. Breadth', '5x Cov. Breadth', '10x Cov. Breadth',
            '20x Cov. Breadth', 'Gini Coefficient']) + '\n'
        )
        for line in summary:
            f.write(sep.join(['{:.8}'.format(i) for i in line]) + '\n')


def plot_lorenz(fig, ax, out_dir, out_name):
    ax.plot([0, 1], [0, 1], '--', color='k', label='Perfect Uniformity')
    ax.set_xlabel('Cumulative Fraction of Genome')
    ax.set_ylabel('Cumulative Fraction of BP Coverage')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_title(out_name)
    ax.legend()
    ax.grid(True)
    
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9)

    out_file = f'{out_name}.lorenz.pdf'
    out_full = os.path.join(out_dir, out_file)
    fig.savefig(out_full, dpi=300)
    


def plot_read_depth_dist(fig, ax, out_dir, out_name):
    ax.set_xlabel('Coverage Depth')
    ax.set_ylabel('Fraction of Coverage $\geq$ depth')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 1])
    # ax.set_xscale('log')
    ax.set_title(out_name)
    ax.legend()
    ax.grid(True)
    
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9)

    out_file = f'{out_name}.readDepthDist.pdf'
    out_full = os.path.join(out_dir, out_file)
    fig.savefig(out_full, dpi=300)


if __name__ == '__main__':
    args = parse_args()
    main(args)