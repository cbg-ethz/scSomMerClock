#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



TICK_FONTSIZE = 16
LABEL_FONTSIZE = 16
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


def plot_wmax_pval(in_dir, out_file=''):
    df = pd.DataFrame(columns=['ADO rate', 'wMax', 'False pos. [%]'])
    for res_dir in os.listdir(in_dir):
        if not res_dir.startswith('res_clock0') or 'bulk' in res_dir:
            continue
        ADO = float(re.search('WGA(0[\.\d]*)', res_dir).group(1))
        filter_dir = os.path.join(in_dir, res_dir, 'minDP5-minGQ1')
        for pTree_dir in os.listdir(filter_dir):
            if not pTree_dir.startswith('poissonTree'):
                continue
            wMax = int(re.search('wMax(\d+)', pTree_dir).group(1))
            wMax_file = os.path.join(filter_dir, pTree_dir, 'poissonTree.summary.tsv')
            with open(wMax_file, 'r') as f:
                res_raw = f.read().strip().split('\n')
            sig_raw =  res_raw[-1].split('\t')[-1].split('/')
            sig = 1 - float(sig_raw[0]) / float(sig_raw[1])

            if np.isclose(wMax, 1000 * ADO, atol=5):
                df.loc[df.shape[0]] = [ADO, r'$(\gamma + \beta)$ 1000', sig]
            else:
                df.loc[df.shape[0]] = [ADO, wMax, sig]

    fig, ax = plt.subplots(figsize=(12, 9))
    sns.lineplot(data=df, x='ADO rate', y='False pos. [%]', hue='wMax',
        palette={101: '#FF7F00', r'$(\gamma + \beta)$ 1000': '#377DB8',
            999: '#E41A1A'}, ax=ax, lw=2
    )

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, default='.',
        help='Base dir for results.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    plot_wmax_pval(args.input, args.output)