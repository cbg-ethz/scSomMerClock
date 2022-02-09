#!/usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
import re

from defaults import *


def generate_ADO_mut_plot(ADO_file, data_file, out_dir):
    ADO = pd.read_csv(ADO_file, sep='\t', index_col=0, encoding='utf-8')
    ADO.drop_duplicates(keep='last', inplace=True)
    if ADO.index.dtype == str:
        ADO.index = [int(i[-4:]) if len(i) > 3 else np.nan \
            for i in ADO.index.values]
    ADO = ADO[~ADO.index.duplicated()].iloc[:,:-1].mean(axis=1)

    ADO_match = re.search('WGA(0[\.\d]*),?(0\.\d*)?-', ADO_file)
    title_str = 'Simulation: '
    try:
        title_str += f'mean={ADO_match.group(1)}, var={ADO_match.group(2)}'
    except IndexError:
        title_str += f'mean={ADO_match.group(1)}'
    except IndexError:
        title_str += '?'
    title_str += '\n'

    data = pd.read_csv(data_file, sep='\t', index_col=0)
    data = data[~data.index.duplicated()]

    plot_data = pd.concat([ADO, data['muts_cells']], axis=1).dropna().values.T

    fig, ax = plt.subplots(figsize=(16, 12))
    plt.plot(plot_data[0], plot_data[1], marker='x', lw=0)

    ax.set_ylabel('# Mutations')
    ax.set_xlabel(f'mean ADO rate')
    ax.tick_params(axis='both', which='major')
    ax.tick_params(axis='both', which='minor')
    plt.title(title_str)

    out_file = os.path.join(out_dir, f'ADO_mut_scatter.png')
    fig.savefig(out_file, dpi=300)
    plt.close()


def generate_est_scatter_plot(est_file, ADO_file, pVal_file, out_dir):
    pVal_df = pd.read_csv(pVal_file, sep='\t', index_col=0)
    est_df = pd.read_csv(est_file, sep='\t', index_col=0)

    ADO_df = pd.DataFrame(
        pd.read_csv(ADO_file, sep='\t', index_col=0).mean(axis=1),
        columns=['FN_true'])

    df_scite = pd.concat(
        [ADO_df, pVal_df['p-value_poissonTree_1.scite'], est_df['Scite_FN'] * 2],
        axis=1) \
        .dropna() \
        .rename({'p-value_poissonTree_1.scite': 'P-value',
            'Scite_FN': 'FN_estimated'}, axis=1)

    df_cellphy = pd.concat(
        [ADO_df, pVal_df['p-value_poissonTree_1.cellphy'], est_df['CellPhy_FN']],
        axis=1) \
        .dropna() \
        .rename({'p-value_poissonTree_1.cellphy': 'P-value',
            'CellPhy_FN': 'FN_estimated'}, axis=1)

    plot_ADO_est(df_scite, os.path.join(out_dir, 'FN_estimate_scite.png'))
    plot_ADO_est(df_cellphy, os.path.join(out_dir, 'FN_estimate_cellphy.png'))


def plot_ADO_est(df, out_file):
    fig, ax = plt.subplots(figsize=(16, 12))

    sc = sns.scatterplot(data=df, x='FN_true', y='FN_estimated', hue='P-value',
        legend=False, ax=ax)


    norm = plt.Normalize(df['P-value'].min(), df['P-value'].max())
    cmap = sns.cubehelix_palette(light=1, as_cmap=True)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cax = fig.add_axes([ax.get_position().x1 - 0.05, 0.15, 0.03, 0.7])
    cbar = ax.figure.colorbar(sm, cax=cax)
    cbar.ax.set_ylabel('P-value', rotation=90)

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.9, hspace=0.33)

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def _get_title(out_file):
    title = out_file.split('_')[-1].replace('.png', '') + ': '
    if 'clock0_' in out_file:
        title += 'Strict Clock '
    else:
        title += 'No Clock '

    try:
        wga_error = re.search('WGA(0\.?\d*)', out_file).group(1)
        if float(wga_error) > 0.1:
            title += '- High Errors '
        elif float(wga_error) > 0:
            title += '- Low Errors'
        elif float(wga_error) == 0:
            title += '- No Errors'
    except AttributeError:
        pass

    try:
        filter_DP = re.search('minDP(\d+)', out_file).group(1)
        filter_GQ = re.search('minGQ(\d+)', out_file).group(1)
        title += f' - Filters (DP>={filter_DP}, GQ>={filter_GQ})'
    except AttributeError:
        pass

    return title


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inDir', help='Input directory')
    parser.add_argument('-a', '--ADO',
        help='ADO file (if not Input directory given.')
    parser.add_argument('-d', '--data',
        help='Summary file (if not Input directory given).')
    parser.add_argument('-e', '--estimates',
        help='Error estimates file (if not Input directory given).')
    parser.add_argument('-p', '--pvalues',
        help='P-values file (if not Input directory given).')
    parser.add_argument('-o', '--outDir', default='', help='Output directory')
    parser.add_argument('-c', '--compress', action='store_true',
        help='Compress output directory')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.inDir:
        ADO_file = os.path.join(args.inDir, 'ADO_overview.tsv')
        data_file = os.path.join(args.inDir, 'data_overview.tsv')
        est_file = os.path.join(args.inDir, 'ADO_estimates.tsv')
        pVal_file = os.path.join(args.inDir, 'final_summary.tsv')
    else:
        ADO_file = args.ADO
        data_file = args.data
        est_file = args.estimates
        pVal_file = args.pvalues
    if not args.outDir:
        if args.inDir:
            args.outDir = os.path.join(args.inDir, 'ADO_plots')
        else:
            args.outDir = os.path.join(os.path.dirname(data_file), 'ADO_plots')
    os.makedirs(args.outDir, exist_ok=True)

    if os.path.exists(est_file):
        generate_est_scatter_plot(est_file, ADO_file, pVal_file, args.outDir)

    generate_ADO_mut_plot(ADO_file, data_file, args.outDir)


    if args.compress:
        import shutil
        shutil.make_archive(args.outDir, 'zip', args.outDir)


