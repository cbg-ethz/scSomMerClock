#!/usr/bin/env python3

import argparse

from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

from defaults import *

NCOL = 3


def plot_legends(args):
    plot_test_legend(args)
    if args.wMax:
        plot_wMax_legend(args)
    if args.power:
        plot_power_legend(args)


def plot_test_legend(args):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.grid(False)
    ax.axis('off')

    handles = []
    labels = []
    tests_all = {
        'general': ['PAUP*', 'poissonDisp'],
        'single-cell': ['cellcoal', 'cellphy', 'scite'],
        'bulk': ['neutrality']
    }
    for test_type, tests in tests_all.items():
        labels.append(fr'$\bf{{{test_type}}}$')
        handles.append(mpatches.Patch(color='white'))

        for test in tests:
            if test not in args.method:
                continue
            labels.append(METHODS[test])
            handles.append(
                mpatches.Patch(color=COLORS[test], label=METHODS[test]))
    ax.legend(handles, labels, ncol=1, frameon=True,
        title=r'$\bf{Tests}$' + '\n' + r'$\bf{Subsample\ size}$' + '\n' + r'$\bf{w_{max}}$')

    labels.append('')
    handles.append(mpatches.Patch(color='white'))
    labels.append(r'$w_{max}$')
    handles.append(mpatches.Patch(color='white'))

    fig.tight_layout()
    if args.output:
        fig.savefig(args.output + '_tests.png', dpi=DPI)
    else:
        plt.show()
    plt.close()


def plot_wMax_legend(args):
    fig, ax = plt.subplots(figsize=(5, 3))

    wMax_vals = [1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    cmap =cm.get_cmap('viridis_r')
    norm = Normalize(vmin=wMax_vals[0], vmax=wMax_vals[-1])

    # handles = []
    # labels = []
    # for wMax in wMax_vals:
    #     if wMax > 300 and wMax < 1000:
    #         handles.append(Line2D([0], [0],
    #             marker='o', color=cmap(norm(wMax)), markersize=5, lw=0))
    #         labels.append('...')
    #     else:
    #         handles.append(Line2D([0], [0],
    #             marker='o', color=cmap(norm(wMax)), markersize=5, lw=0))
    #         labels.append(wMax)
    # ax.legend(handles, labels, ncol=2, frameon=True, title=r'$\bf{w_{max}}$')

    cb = ColorbarBase(
        ax,
        orientation='horizontal',
        cmap=cmap,
        norm=norm,
        ticks=[1, 500, 1000],

    )
    cb.ax.tick_params(labelsize=18)

    fig.tight_layout()
    if args.output:
        fig.savefig(args.output + '_wMax.png', dpi=DPI)
    else:
        plt.show()
    plt.close()


def plot_power_legend(args):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.grid(False)
    ax.axis('off')

    handles = []
    labels = []
    for ss in args.subsamples:
        labels.append(f'{ss} / {args.total_cells}')
        handles.append(Line2D(
            [0], [0], color='grey', markersize=0, lw=2, linestyle=LINE_STYLE[ss]))
    ax.legend(handles, labels, ncol=1, frameon=True, title=r'$\bf{subsamples}$')

    fig.tight_layout()
    if args.output:
        fig.savefig(args.output + '_power.png', dpi=DPI)
    else:
        plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-w', '--wMax', action='store_true',
        help='Plot w_max legend.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellcoal', 'cellphy', 'scite', 'poissonDisp', 'PAUP*', 'neutrality', 'mobster'],
        default=['cellcoal', 'cellphy', 'poissonDisp', 'PAUP*', 'neutrality', 'mobster'],
        help='Plot method legend.')
    parser.add_argument('-p', '--power', action='store_true',
        help='Plot power legend.')
    parser.add_argument('-ss', '--subsamples', nargs='+', type=int,
        choices=[10, 30, 50, 70, 90], default=[10, 50, 90],
        help='# Cells subsampled. Default = [10, 50, 90].')
    parser.add_argument('-t', '--total_cells', type=int, default=100,
        help='Total number of simulated cells. Default = 100.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    plot_legends(args)