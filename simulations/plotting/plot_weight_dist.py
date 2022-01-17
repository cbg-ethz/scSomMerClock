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
from matplotlib.patches import Patch
from matplotlib.transforms import ScaledTranslation
from scipy.stats import chi2
from matplotlib.colors import to_rgba


COLORS = [
     # Blue     # Green    # Red      # Orange   # Purple
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', # light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]

TICK_FONTSIZE = 16
LABEL_FONTSIZE = 16

vis_names = {'cellcoal': 'True', 'scite': 'SCITE', 'cellphy': 'CellPhy'}
colors = {
    'True': (0.945, 0.4, 0.627),
    'CellPhy': (0.443, 0.396, 0.776),
    'SCITE': (0.02, 0.604, 0.173)
}
colors = {101: '#FF7F00', r'$(\gamma + \beta)$ 1000': '#377DB8', 999: '#E41A1A'}

def generate_weights_plot(in_files, out_file=None):
    data = {}
    legend_elements = []
    for i, in_file in enumerate(in_files):
        df_in = pd.read_csv(in_file, sep='\t', engine='python', index_col='run',
        skipfooter=1)
        weights_raw = df_in.values[:, -1]
        weights = []
        for run_weights in weights_raw:
            weights.append(np.array(run_weights.split(','), dtype=float))

        _, wMax_raw, tree, = in_file.split(os.sep)[-2].split('_')
        wMax = int(wMax_raw.replace('wMax', ''))

        ADO = float(re.search('WGA(0[\.\d]*)', in_file).group(1))
        if np.isclose(wMax, 1000 * ADO, atol=5):
            label = r'$(\gamma + \beta)$ 1000'
        else:
            label = wMax

        data[label] = np.array(weights).flatten()

        legend_elements.append(Patch(fc=to_rgba(colors[label], 0.5),
            ec=colors[label], label=label, linewidth=1))

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

    # Same as p-value plot for easy manual merging
    fig, ax = plt.subplots(figsize=(15, 9))

    sns.histplot(data,
        element='poly', stat='density', kde=False,
        common_norm=False, fill=True,
        bins=100, multiple='dodge',
        kde_kws={'cut': 0, 'clip': (0, 1)},
        line_kws={'lw': 3, 'alpha': 0.75},
        alpha=0.4, palette=colors,
        legend=False, ax=ax
    )
    ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax.set_xlabel('Poisson Tree weights')

    plt.legend(handles=legend_elements, title='Tree:', bbox_to_anchor=(1.1, 1),
        loc='upper left', borderaxespad=0.)

    fig.subplots_adjust(left=0.11, bottom=0.15, right=0.75, top=0.95)

    if not out_file:
        out_dir = os.sep.join(in_file.split(os.sep)[:-2])
        out_file = os.path.join(out_dir, 'PoissonTree_weights.pdf')
    fig.savefig(out_file, dpi=300)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input directory.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    in_files = []
    for res_dir in os.listdir(args.input):
        if res_dir.startswith('poissonTree'):
            in_files.append(
                os.path.join(args.input, res_dir, 'poissonTree.summary.tsv'))
    generate_weights_plot(in_files, args.output)
