#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from scipy.stats import chi2


COLORS = [
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', #light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]

TICK_FONTSIZE = 12
LABEL_FONTSIZE = 16

vis_names = {
    'poisson': 'Poisson',
    'paup': 'PAUP*',
    'poissontree': 'Poisson + Tree',
    'cellcoal': 'True',
    '-': '-',
    'scite': 'SCITE',
    'cellphy': 'CellPhy',

}


def generate_pval_plot(in_file, out_file=None, bins=0):
    df_in = pd.read_csv(in_file, sep='\t', engine='python', index_col='run',
        skipfooter=1)
    n = df_in.shape[0]

    if not bins:
        bins = 100 # int(np.sqrt(n))

    df = pd.DataFrame(columns=['P-value', 'Method', 'Tree'])
    df_lambda = pd.DataFrame(columns=['Lambda', 'Tree'])
    for rel_col in [i for i in df_in.columns if 'p-value' in i]:
        df_new = pd.DataFrame(df_in[rel_col]).rename({rel_col: 'P-value'}, axis=1)
        try:
            method, tree = rel_col.lower().replace('p-value_', '').split('.')
        except ValueError:
            method = rel_col.lower().replace('p-value_', '')
            tree = '-'

        df_new['Method'] = vis_names[method]
        df_new['Tree'] = vis_names[tree]
        df = df.append(df_new, ignore_index=True)

        if vis_names[method] == 'Poisson + Tree':
            lambda_col = rel_col.replace('p-value_', '-2logLR_')
            df_lambda_new = pd.DataFrame(
                df_in[lambda_col]).rename({lambda_col: 'Lambda'}, axis=1)
            df_lambda_new['Tree'] = vis_names[tree]
            df_lambda = df_lambda.append(df_lambda_new, ignore_index=True)

    sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
    sns.set_context('paper',
        rc={'font.size': TICK_FONTSIZE,
            'axes.titlesize': LABEL_FONTSIZE,
            'axes.labelsize': LABEL_FONTSIZE,
            'axes.titlesize': LABEL_FONTSIZE,
            'axes.labelticksize': LABEL_FONTSIZE,
            'lines.linewidth': 2,
            'legend.fontsize': LABEL_FONTSIZE,
            'legend.title_fontsize':  LABEL_FONTSIZE,
            'xtick.major.size':  TICK_FONTSIZE*2,
            'ytick.major.size':  TICK_FONTSIZE,
        })

    colors = {
        'SCITE': '#F6577C',
        'CellPhy': '#47AF09',
        'True': '#153386',
        '-': '#C88B0A'
    }
    plot_pvals(df, out_file, bins, colors)
    plot_test_statistic(df_lambda, out_file, bins, colors)


def plot_pvals(df, out_file, bins, colors):

    row_order = ['Poisson', 'PAUP*', 'Poisson + Tree']
    dp = sns.displot(df, x='P-value', hue='Tree', row='Method',
        element='bars', stat='density', kde=False,
        common_norm=False, fill=True, bins=bins,
        kde_kws={'cut': 0, 'clip': (0, 1)},
        line_kws={'lw': 3, 'alpha': 0.75},
        palette=colors, lw=2, alpha=0.5,
        height=4, aspect=3,
        hue_order=['-', 'SCITE', 'True', 'CellPhy'],
        row_order=row_order,
        facet_kws={'sharey': False})

    ideal_y = 1
    for ax_no, method in enumerate(row_order):
        y_lim = max([i.get_height() for i in dp.axes[ax_no][0].patches]) * 1.05
        # y_lim = min(1, ideal_y * 5)
        dp.axes[ax_no][0].axhline(ideal_y, ls='--', c='red', lw=2)

        dp.axes[ax_no][0].set_ylim([0, y_lim])
        dp.axes[ax_no][0].set_xlim([0, 1])

        dp.axes[ax_no][0].tick_params(axis='both', which='major',
            labelsize=TICK_FONTSIZE)

    # Rugplot for 1 row: poisson
    height = -0.03
    poisson_data = df[(df['Tree'] == '-') & (df['Method'] == 'Poisson')]
    sns.rugplot(data=poisson_data, x='P-value', height=height, clip_on=False,
        color=colors['-'], ax=dp.axes[0][0], linewidth=0.2)
    dp.axes[0][0].add_collection(
        LineCollection(np.array([[[0.05, 0],[0.05, height]]]),
        transform=dp.axes[0][0].get_xaxis_transform(),
        clip_on=False, color='red', ls='-', linewidth=2)
    )
    sgf_perc = np.percentile(poisson_data['P-value'].values, 5)
    dp.axes[0][0].axvline(sgf_perc, ls='--', color=colors['-'], lw=2)
    dp.axes[0][0].text(sgf_perc + 0.005, dp.axes[0][0].get_ylim()[1] / 2,
        f'5th percentile ({sgf_perc:.3f})', color=colors['-'])
    # dp.axes[0][0].add_collection(
    #     LineCollection(np.array([[[sgf_perc, 0], [sgf_perc, height]]]),
    #     transform=dp.axes[0][0].get_xaxis_transform(),
    #     clip_on=False, alpha=1, color='black', ls='-', linewidth=2)
    # )

    # Rugplot for 2 and 3 row: paup and poisson + tree
    for ax_no, method in [(1, 'PAUP*'), (2, 'Poisson + Tree')]:
        for i, tree in enumerate(['True', 'SCITE', 'CellPhy']):
            data = df[(df['Tree'] == tree) & (df['Method'] == method)] \
                ['P-value'].values
            if data.size == 0:
                continue

            segs = np.stack((np.c_[data, data],
                np.c_[np.zeros_like(data) + height * i,
                        np.zeros_like(data) + height * (i + 1)]),
                    axis=-1)
            lc = LineCollection(segs,
                transform=dp.axes[ax_no][0].get_xaxis_transform(),
                clip_on=False, color=colors[tree], linewidth=0.025, alpha=0.75)
            dp.axes[ax_no][0].add_collection(lc)

            sgf_perc = np.percentile(data, 5)
            dp.axes[ax_no][0].axvline(sgf_perc, ls='--', color=colors[tree], lw=2)
            dp.axes[ax_no][0].text(sgf_perc + 0.01,
                dp.axes[ax_no][0].get_ylim()[1] / 6 * (i + 3),
                f'5th percentile ({sgf_perc:.3f})', color=colors[tree])

            # sgf_line = np.array([[[sgf_perc, height * i],
            #     [sgf_perc, height * (i + 1)]]])
            # dp.axes[ax_no][0].add_collection(
            #     LineCollection(sgf_line,
            #     transform=dp.axes[ax_no][0].get_xaxis_transform(),
            #     clip_on=False, alpha=1, color='black', ls='-', linewidth=2)
            # )

        pval_line = np.array([[[0.05, 0],[0.05, height * 3]]])
        dp.axes[ax_no][0].add_collection(
            LineCollection(pval_line,
            transform=dp.axes[ax_no][0].get_xaxis_transform(),
            clip_on=False, color='red', ls='-', linewidth=2)
        )


    dp.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.85, top=0.9, hspace=0.33)
    dp.fig.suptitle(_get_title(in_file), fontsize=LABEL_FONTSIZE*1.25)

    if out_file:
        dp.fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_test_statistic(df, out_file, bins, colors):
    dof = 29 # TODO <NB> remove hardcoding

    if not bins:
        bins = 100

    row_order = ['SCITE', 'True', 'CellPhy']
    # color = tuple([colors[i] for i in row_order])
    # color = None

    dp = sns.displot(df, x='Lambda', row='Tree',
        element='bars', stat='density', kde=False,
        common_norm=False, common_bins=False, fill=True, bins=bins,
        kde_kws={'cut': 0, 'clip': (0, 1)},
        line_kws={'lw': 3, 'alpha': 0.75},
        lw=2, alpha=0.5,
        height=3, aspect=3,
        row_order=['SCITE', 'True', 'CellPhy'])

    x = np.linspace(chi2.ppf(0.001, dof), chi2.ppf(0.999, dof), bins)
    y = chi2.pdf(x, dof)
    for ax_no, method in enumerate(row_order):
        dp.axes[ax_no][0].plot(x, y, ls='--', color='r', lw=2,
                label=f'Chi2({dof}) pdf')
        for child in dp.axes[ax_no][0].properties()['children'][:-1]:
            if child.__str__().startswith('Rectangle'):
                child.set_color(colors[method])

    if out_file:
        dp.fig.savefig(out_file.replace('.pdf', '.Lambda.pdf'), dpi=300)
    else:
        plt.show()
    plt.close()



def _get_title(in_file):
    title = ''
    if 'clock0_' in in_file:
        title += 'Strict Clock '
    else:
        title += 'No Clock '

    try:
        wga_error = re.search('WGA(0\.?\d*)', in_file).group(1)
        if float(wga_error) > 0.1:
            title += '- High Errors '
        elif float(wga_error) > 0:
            title += '- Low Errors'
        elif float(wga_error) == 0:
            title += '- No Errors'
    except AttributeError:
        pass

    try:
        filter_DP = re.search('minDP(\d+)', in_file).group(1)
        filter_GQ = re.search('minGQ(\d+)', in_file).group(1)
        filter_singleton = 'noSingletons' in in_file
        if int(filter_DP) == 0 and int(filter_GQ) == 0 and not filter_singleton:
            title += ' - No Filters'
        else:
            title += f' - Filters (DP>={filter_DP}, GQ>={filter_GQ}'
            if filter_singleton:
                title += ', no Singletons)'
            else:
                title += ')'
    except AttributeError:
        pass

    return title


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-b', '--bins', type=int, default=0,
        help='Bin number in histogram. Default = sqrt(n).')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if os.path.isdir(args.input):
        rel_files = [os.path.join(args.input, i) for i in os.listdir(args.input) \
            if i.startswith('final_summary') and i.endswith('.tsv')]
        for in_file in sorted(rel_files):
            out_file = f'{in_file}.pdf'
            print(f'Processing: {in_file}')
            generate_pval_plot(in_file, out_file, args.bins)
    else:
        if not args.output:
            args.output = f'{args.input}.pdf'
        generate_pval_plot(args.input, args.output, args.bins)