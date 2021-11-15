#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re


TICK_FONTSIZE = 16
LABEL_FONTSIZE = 20

COLORS = [
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', #light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]


def plot_ADO_distr(ADO, out_file, title_str=''):
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.hist(ADO.values.flatten(), density=True, histtype='stepfilled',
        color=COLORS[1], label='per cell', alpha=0.9)
    ax.hist(ADO.mean(axis=1).values, density=True, histtype='stepfilled',
        color=COLORS[0], label='per run', alpha=0.9)
    ax.set_xlabel('ADO', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Density', fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
    plt.legend(fontsize=TICK_FONTSIZE)
    plt.title(title_str, fontsize=LABEL_FONTSIZE)

    fig.savefig(out_file, dpi=300)
    plt.close()


def generate_pval_plot(ADO_file, data_file, out_dir):
    ADO = pd.read_csv(ADO_file, sep='\t', index_col=0, encoding='utf-8')
    ADO.drop_duplicates(keep='last', inplace=True)
    if ADO.index.dtype == str:
        ADO.index = [int(i[-4:]) if len(i) > 3 else np.nan \
            for i in ADO.index.values]
    ADO = ADO[~ADO.index.duplicated()]
    #import pdb; pdb.set_trace()

    ADO_match = re.search('WGA(0[\.\d]*),?(0\.\d*)?-', ADO_file)
    title_str = 'Simulation: '
    try:
        title_str += f'mean={ADO_match.group(1)}, var={ADO_match.group(2)}'
    except IndexError:
        title_str += f'mean={ADO_match.group(1)}'
    except IndexError:
        title_str += '?'
    title_str += '\n'

    plot_ADO_distr(ADO, os.path.join(out_dir, 'ADO_distr.png'), title_str)

    data = pd.read_csv(data_file, sep='\t', index_col=0)

    methods = [i for i in data.columns if i.startswith('p-value')]

    y_all = [('Mean', ADO.mean(axis=1)), ('StD', ADO.std(axis=1))]
    for method in methods:
        method_str = method.replace('p-value', '').strip('_')
        x = data[method]
        for y_label, y in y_all:
            try:
                plot_data = pd.concat([y, x], axis=1).dropna().values.T
            except:
                import pdb; pdb.set_trace()

            fig, ax = plt.subplots(figsize=(16, 12))
            plt.plot(plot_data[0], plot_data[1], marker='x', lw=0)

            ax.set_ylabel('p-Value', fontsize=LABEL_FONTSIZE)
            ax.set_xlabel(f'ADO {y_label}', fontsize=LABEL_FONTSIZE)
            ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
            ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
            plt.title(title_str + f'Method: {method_str}', fontsize=LABEL_FONTSIZE)

            out_file = os.path.join(out_dir, f'ADO_pVal_{method_str}_{y_label}.png')
            fig.savefig(out_file, dpi=300)
            plt.close()

        # 3D plot
        d_data = pd.concat([y_all[0][1], y_all[1][1], x], axis=1).dropna().values.T
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(projection='3d')
        ax.scatter(d_data[0], d_data[1], d_data[2], marker='x')
        ax.set_xlabel('ADO Mean')
        ax.set_ylabel('ADO StD')
        ax.set_zlabel('p-Value')
        ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
        plt.title(title_str + f'Method: {method_str}', fontsize=LABEL_FONTSIZE)

        out_file = os.path.join(out_dir, f'ADO_pVal_{method_str}_3D.png')
        fig.savefig(out_file, dpi=300)
        plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inDir', help='Input directory')
    parser.add_argument('-a', '--ADO',
        help='ADO file (if not Input directory given.')
    parser.add_argument('-s', '--summary',
        help='Summary file (if not Input directory given.')
    parser.add_argument('-o', '--outDir', default='', help='Output directory')
    parser.add_argument('-c', '--compress', action='store_true',
        help='Compress output directory')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    if args.inDir:
        ADO_file = os.path.join(args.inDir, 'ADO_overview.tsv')
        data_file = os.path.join(args.inDir, 'final_summary.tsv')
    else:
        ADO_file = args.ADO
        data_file = args.summary
    if not args.outDir:
        if args.inDir:
            args.outDir = os.path.join(args.inDir, 'ADO_plots')
        else:
            args.outDir = os.path.join(os.path.dirname(data_file), 'ADO_plots')
        os.makedirs(args.outDir, exist_ok=True)

    generate_pval_plot(ADO_file, data_file, args.outDir)
    if args.compress:
        import shutil
        shutil.make_archive(args.outDir, 'zip', args.outDir)


