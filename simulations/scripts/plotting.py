#!/usr/bin/env python3

import argparse
from scipy.stats import linregress, gamma
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter


TICK_FONTSIZE = 16
LABEL_FONTSIZE = 20


def generate_mrbayes_plots(in_file, out_file, regress=False):
    def get_ratio(ev_str):
        success, total = ev_str.split('/')
        return float(success) / float(total)

    df = pd.read_csv(in_file, sep='\t')
    df['ratio'] = df['Evidence'].apply(lambda x: get_ratio(x))
    df['runtime'] = df['Avg. runtime [secs]'].apply(lambda x: x / 60 / 60)

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 1)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[2, 0])

    ax0.errorbar(df['steps'], df['Avg. H_0:clock'], yerr= df['Std. H_0:clock'],
        label=r'$H_0$', color='#1F78B4', capsize=4, ls='--', marker='x')
    ax0.errorbar(df['steps'], df['Avg. H_1:noClock'], yerr= df['Std. H_1:noClock'],
        label=r'$H_1$', color='#33A02C', capsize=4, ls='--', marker='x')
    ax0.set_ylabel('Log-likelihood', fontsize=LABEL_FONTSIZE)
    ax0.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
    ax0.legend(fontsize=TICK_FONTSIZE)

    ax1.errorbar(df['steps'], df['Avg. 2log_e(B_01)'],
        yerr= df['Std. 2log_e(B_01)'], color='#FF7F00', capsize=4, ls='--',
        marker='x')
    ax1.set_ylabel(r'$2 log_e(B_{01})$', fontsize=LABEL_FONTSIZE)
    ax1.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))

    ax2.plot(df['steps'], df['ratio'], color='#E31A1C', ls='', marker='x',
        label='real')

    if regress:
        reg_line = linregress(df['steps'], df['ratio'])
        reg_x_alpha = (0.05 - reg_line.intercept) / reg_line.slope
        reg_x = np.linspace(df['steps'].min(), np.ceil(reg_x_alpha / 1e7) * 1e7, 20)
        reg_y = reg_line.intercept + reg_line.slope * reg_x

        ax2.plot(reg_x, reg_y, color='#6A3D9A', ls='--', marker='x',
            label=f'{reg_line.intercept:.2f} + {reg_line.slope:.2E} * x    '
                f'($R^2$={reg_line.rvalue:.2f})')
        ax2.axhline(0.05, ls=':')
        ax2.axvline(reg_x_alpha, ls=':')
        ax2.text(reg_x_alpha, 0.05, f'({reg_x_alpha:.2E}, 0.05)', va='top',
            ha='right', fontsize=TICK_FONTSIZE)
        ax2.legend(fontsize=TICK_FONTSIZE)

    ax2.set_ylabel(r'$H_0$ rejected [%]', fontsize=LABEL_FONTSIZE)
    ax2.set_xlabel('MCMC steps', fontsize=LABEL_FONTSIZE)
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
    ax2.set_ylim(-.05, 1.05)

    reg_line2 = linregress(df['steps'], df['runtime'])
    def mcmc_to_runtime(x):
        return reg_line2.intercept + reg_line2.slope * x
    
    secax = ax2.secondary_xaxis('top', functions=(mcmc_to_runtime, mcmc_to_runtime))
    secax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    secax.set_xlabel('Runtime [h]')   

    fig.suptitle('Simulations: no WGS, no NGS', fontsize=LABEL_FONTSIZE * 1.5)

    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92,
        hspace=0.5)
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def generate_pval_plot(in_file, out_file, p_val_filter=1):
    df = pd.read_csv(in_file, sep='\t', engine='python', skipfooter=1)
    fig, ax = plt.subplots(figsize=(16, 12))

    p_vals = df[df['p-value'] <= p_val_filter]['p-value']
    if p_val_filter < 1:
        print('Removed p-values: {} / {}' \
            .format(df.shape[0] - p_vals.size, df.shape[0]))
    # p_vals = df['p-value']

    ax.hist(p_vals, bins=100, range=(0, 1))
    ax.set_ylabel(f'counts (n={p_vals.size})', fontsize=LABEL_FONTSIZE)
    ax.set_xlabel('p-values', fontsize=LABEL_FONTSIZE)
    ax.set_xlim([-0.01, 1.01])
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92,
        hspace=0.5)
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def generate_gamma_plot(alphas, out_file):
    fig, ax = plt.subplots(figsize=(16, 12))

    y_max_t = 0
    x_max_t = 0
    for alpha in sorted(alphas):
        scale = 1 / alpha
        x_min = gamma.ppf(0.001, alpha, scale=scale)
        x_max = gamma.ppf(0.999, alpha, scale=scale)
        x = np.linspace(x_min, x_max, 1000)
        y = gamma.pdf(x, alpha, scale=scale)
        ax.plot(x, y, '-', label=fr'$\alpha$ = {alpha}')
        y_max_t = max(y_max_t, y.max())
        x_max_t = max(x_max_t, x.max())

    ax.axvline(1, ls='--', c='black')

    ax.set_ylabel(r'Gamma(x; $\alpha$, $1 / \alpha$)', fontsize=LABEL_FONTSIZE)
    ax.set_ylim([0, y_max_t])
    ax.set_xlabel('x', fontsize=LABEL_FONTSIZE)
    ax.set_xlim([0, x_max_t])

    ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=TICK_FONTSIZE)
    ax.set_title('PDF', fontsize=LABEL_FONTSIZE * 1.5)

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input file')
    parser.add_argument('-f', '--format', type=str, default='pval',
        choices=['pval', 'mrbayes', 'gamma'], help='What type of plot to generate')
    parser.add_argument('-a', '--alpha', type=float, default=[1], nargs='+',
        help='Alpha parameters of Gamma distribution with mean 1. Default = 1')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
   
    if args.format == 'mrbayes':
        generate_mrbayes_plots(args.input,args.output)
    elif args.format == 'pval':
        generate_pval_plot(args.input, args.output)
    elif args.format == 'gamma':
        generate_gamma_plot(args.alpha, args.output)
    else:
        raise IOError('Unknown format type: {}'.format(args.format))