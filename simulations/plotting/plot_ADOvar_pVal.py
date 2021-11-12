#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


TICK_FONTSIZE = 16
LABEL_FONTSIZE = 20


def plot_ADO_dist(ADO, out_file):
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.hist(ADO.values.flatten(), density=True, histtype='stepfilled')
    ax.set_xlabel('ADO (per cell)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Density', fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
    fig.savefig(out_file, dpi=300)
    plt.close()


def generate_pval_plot(ADO_file, data_file, out_dir):
    ADO = pd.read_csv(ADO_file, sep='\t', index_col=0)
    ADO.drop_duplicates(keep='last', inplace=True)

    plot_ADO_dist(ADO, os.path.join(out_dir, 'ADO_dist.png'))

    data = pd.read_csv(data_file, sep='\t', index_col=0)

    methods = [i for i in data.columns if i.startswith('p-value')]

    y_all = [('Mean', ADO.mean(axis=1)), ('StD', ADO.std(axis=1))]
    for method in methods:
        method_str = method.replace('p-value', '').strip('_')
        x = data[method]
        for y_label, y in y_all:
            plot_data = pd.concat([y, x] ,axis=1).dropna().values.T

            fig, ax = plt.subplots(figsize=(16, 12))
            plt.plot(plot_data[0], plot_data[1], marker='x', lw=0)

            ax.set_ylabel('p-Value', fontsize=LABEL_FONTSIZE)
            ax.set_xlabel(f'ADO {y_label}', fontsize=LABEL_FONTSIZE)
            ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
            ax.tick_params(axis='both', which='minor', labelsize=TICK_FONTSIZE)
            plt.title(f'Method: {method_str}')

            out_file = os.path.join(out_dir, f'ADO_pVal_{method_str}_{y_label}.png')
            fig.savefig(out_file, dpi=300)
            plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inDir', help='Input directory')
    parser.add_argument('-o', '--outDir', default='', help='Output directory')
    parser.add_argument('-c', '--compress', action='store_true',
        help='Compress output directory')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    ADO_file = os.path.join(args.inDir, 'ADO_overview.tsv')
    data_file = os.path.join(args.inDir, 'final_summary.tsv')
    if not args.outDir:
        args.outDir = os.path.join(args.inDir, 'ADO_plots')
        os.makedirs(args.outDir, exist_ok=True)

    generate_pval_plot(ADO_file, data_file, args.outDir)
    if args.compress:
        import shutil
        shutil.make_archive(args.outDir, 'zip', args.outDir)


