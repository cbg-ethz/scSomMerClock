#!/usr/bin/env python3

import argparse
from scipy.stats import beta
import seaborn as sns
import matplotlib.pyplot as plt


COLORS = [
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', #light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]

TICK_FONTSIZE = 8
LABEL_FONTSIZE = 8



def plot_betabinomial(mean, var, bins=0, density=False):
    if (var >= mean*(1 - mean)):
        raise RuntimeError('Variance has to be smaller than: mean * (1 - mean)')

    n = 1000
    if not bins:
        bins = 100
    a = mean * ((mean * (1 - mean) / var) - 1.0)
    b = (1 - mean) * ((mean * (1 - mean) / var) - 1.0)

    r = beta.rvs(a,b, size=n)

    if density:
        sns.set_style('whitegrid')
        d = sns.kdeplot(r)
        d.set_xlim([0, 1])
    else:
        fig, ax = plt.subplots(1, 1)
        ax.hist(r, density=True, histtype='stepfilled')

    plt.show()
    plt.close()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=float, help='Mean.')
    parser.add_argument('-v', type=float, help='Variance.')
    parser.add_argument('-b', '--bins', type=int, default=0,
        help='Bin number. Default = 100.')
    parser.add_argument('-d', action='store_true',
        help='Plot density instead of histogram.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    plot_betabinomial(args.m, args.v, args.bins, args.d)