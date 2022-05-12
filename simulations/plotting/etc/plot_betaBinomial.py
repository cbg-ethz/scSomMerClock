#!/usr/bin/env python3

import argparse
from scipy.stats import beta

from defaults import *


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

    print(f'Beta: {beta.mean(a, b):.2f} =/- {beta.std(a, b):.3f}')

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