#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
from scipy.stats import chi2


TICK_FONTSIZE = 8
LABEL_FONTSIZE = 8
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'font.size': TICK_FONTSIZE,
        'axes.titlesize': LABEL_FONTSIZE,
        'axes.labelsize': LABEL_FONTSIZE,
        'axes.titlesize': LABEL_FONTSIZE,
        'axes.labelticksize': LABEL_FONTSIZE,
        'lines.linewidth': 1,
        'legend.fontsize': LABEL_FONTSIZE,
        'legend.title_fontsize':  LABEL_FONTSIZE,
        'xtick.major.size':  TICK_FONTSIZE*2,
        'ytick.major.size':  TICK_FONTSIZE,
    })

COLORS = [
     # Blue     # Green    # Red      # Orange   # Purple
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', # light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]
vis_names = {
    'poissondisp': 'Poisson Dispersion',
    'paup': 'PAUP*',
    'poissontree': 'Poisson Tree',
    'cellcoal': 'True Tree',
    '-': '-',
    'scite': 'SCITE Tree',
    'cellphy': 'CellPhy Tree',
}
colors = {
    'Poisson Dispersion':  '#994EA3', # purple
    'PAUP* + True Tree': '#4FAF4A', # green
    'PAUP* + CellPhy Tree': '#177512', # darker green
    'PAUP* + SCITE Tree': '#AAE3A7', #  lighter green
}
poisson_colors = { # red, blue, orange
    'True Tree': ['#E41A1A', '#377DB8', '#FF7F00'], # normal
    'CellPhy Tree': ['#8C0000', '#094D85', '#9B4D00'], # darker
    'SCITE Tree': ['#F04949', '#84B5DE', '#FFB164'] #brigher
}


def generate_pval_plot(in_file, out_file=None, skip_poisson=False):
    df_in = pd.read_csv(in_file, sep='\t', engine='python', index_col='run',
        skipfooter=1)

    df = pd.DataFrame(columns=['P-value', 'Method'])
    df_lambda = pd.DataFrame(columns=['Lambda', 'Tree'])
    for rel_col in sorted([i for i in df_in.columns if 'p-value' in i]):
        df_new = pd.DataFrame(df_in[rel_col]).rename({rel_col: 'P-value'}, axis=1)
        rel_col_short = rel_col.lower().replace('p-value_', '')

        try:
            method, tree = rel_col_short.split('.')
        except ValueError:
            if skip_poisson:
                continue
            method = 'poissondisp'
            method_full = vis_names[method]
        else:
            method = method.split('_')[0]
            method_full = f'{vis_names[method]} + {vis_names[tree]}'
            if method == 'poissontree':
                weights = rel_col_short.split('.')[0].split('_')[1] \
                    .replace('wmax', '')
                method_full += r' + $w_{max}$ = ' + weights
                colors[method_full] = poisson_colors[vis_names[tree]].pop()

        df_new['Method'] = method_full
        df = df.append(df_new, ignore_index=True)

        if vis_names[method] == 'Poisson Tree':
            lambda_col = rel_col.replace('p-value_', '-2logLR_')
            df_lambda_new = pd.DataFrame(
                df_in[lambda_col]).rename({lambda_col: 'Lambda'}, axis=1)
            df_lambda_new['Tree'] = vis_names[tree]
            df_lambda = df_lambda.append(df_lambda_new, ignore_index=True)

    plot_pvals(df, out_file)
    if not df_lambda.empty:
        plot_test_statistic(df_lambda, out_file, colors)


def plot_pvals(df, out_file):
    hue_order = sorted(df.Method.unique())

    fig = plt.figure(figsize=(9, 3))
    gs = GridSpec(1, 4)
    ax1 = fig.add_subplot(gs[0, 0:3])
    ax2 = fig.add_subplot(gs[0, 3])

    dp = sns.histplot(df, x='P-value', hue='Method',
        element='poly', stat='density', kde=False,
        common_norm=False, fill=False,
        binwidth=0.01, binrange=(0, 1), multiple='dodge',
        kde_kws={'cut': 0, 'clip': (0, 1)},
        line_kws={'lw': 3, 'alpha': 0.66},
        palette=colors, alpha=0.5,
        hue_order=hue_order, legend=False,
        ax=ax1,
    )

    # axc1 = dp.axes[0][0]
    # fig = dp.fig

    ax1.set_ylim([min(-4), None])
    ax1.set_xlim([0, 1])
    ax1.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)

    # Add rugplots
    height = 0.03
    legend_elements = []
    h_sum = []
    for i, method in enumerate(hue_order):
        data = df[df['Method'] == method]['P-value'].values
        if data.size == 0:
            continue

        segs = np.stack((np.c_[data, data],
            np.c_[np.zeros_like(data) + 1 + height * i,
                    np.zeros_like(data) + 1 + height * (i + 1)]),
                axis=-1)
        lc = LineCollection(segs,
            transform=ax1.get_xaxis_transform(),
            clip_on=False, color=colors[method], linewidth=0.05, alpha=0.75)
        ax1.add_collection(lc)

        h_sum.append(np.mean(data <= 0.05))
        legend_elements.append(
            Line2D([0], [0], color=colors[method], lw=2, label=method))

    sns.barplot(x=hue_order, y=h_sum, ax=ax2, hue_order=hue_order, palette=colors)
    ax2.set_ylim([0, 1])
    ax2.set_ylabel('P-values $\leq 0.05$ [%]')
    ax2.set_xlabel('Method')
    ax2.axhline(0.05, ls='--', color='red', lw=1)
    ax2.set_xticklabels([])

    plt.legend(handles=legend_elements, title='Method:', bbox_to_anchor=(1.1, 1),
        loc='upper left', borderaxespad=0.)

    sns.despine(bottom=True)
    fig.subplots_adjust(left=0.075, bottom=0.25, right=0.7, top=0.75, wspace=0.75,
        hspace=0)
    fig.suptitle(_get_title(in_file), fontsize=LABEL_FONTSIZE*1.25)

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_test_statistic(df, out_file, colors):
    dof = 29 # TODO <NB> remove hardcoding
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
        ADO_error = re.search('WGA(0\.?\d*)', in_file).group(1)
        title += f'- ADO = {ADO_error} '
    except AttributeError:
        pass

    return title


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-spd', '--skip_poisson', action='store_true',
        help='Skip Poisson distribution test.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if os.path.isdir(args.input):
        rel_files = [os.path.join(args.input, i) for i in os.listdir(args.input) \
            if ('final_summary' in i) and i.endswith('.tsv')]

        for in_file in sorted(rel_files):
            out_file = f'{in_file}.pdf'
            print(f'Processing: {in_file}')
            generate_pval_plot(in_file, out_file, args.skip_poisson)
    else:
        if not args.output:
            args.output = f'{args.input}.pdf'
        generate_pval_plot(args.input, args.output, args.skip_poisson)