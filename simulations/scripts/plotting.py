#!/usr/bin/env python3

import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter


TICK_FONTSIZE = 16
LABEL_FONTSIZE = 20


def generate_mrbayes_plots(in_file, out_file, regress=False):
    from scipy.stats import linregress


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


def simulate_poisson_dispersion(c=30, n=10000):
    from scipy.stats import poisson
    from scipy.stats.distributions import chi2

    # mu = np.clip(np.random.normal(300, 200, size=n), 1e-6, np.inf)
    mu = np.clip(np.random.exponential(100, size=n), 1e-6, np.inf)

    X = poisson.rvs(mu, size=(c, n)).T
    X_mean = X.mean(axis=1, keepdims=True)

    ll_clock = -np.sum(X * np.log(X_mean) - X_mean, axis=1)
    ll_uncon = -np.sum(X * np.log(X) - X, axis=1)
    LR = 2 * (ll_clock - ll_uncon)
    p_vals = chi2.sf(LR, c - 1)

    return p_vals


def _plot_pvals(p_vals, bin_no=100):
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.hist(p_vals, bins=bin_no, range=(0, 1))
    ax.axhline(p_vals.size / bin_no, ls='--', c='red')
    ax.set_ylabel(f'counts (n={p_vals.size})', fontsize=LABEL_FONTSIZE)
    ax.set_xlabel('p-values', fontsize=LABEL_FONTSIZE)
    ax.set_xlim([-0.01, 1.01])
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92, hspace=0.5)
    return fig


def _plot_muts(muts, bin_no=100):
    mu = np.mean(muts)
    half_stdev = np.sqrt(np.var(muts)) / 2

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.hist(muts, bins=bin_no, range=(-1, max(muts) + 10 ))
    ax.axvline(mu, ls='-', c='black')
    ax.axvline(mu - half_stdev, ls='--', c='black')
    ax.axvline(mu + half_stdev, ls='--', c='black')
    ax.set_xlim([-1, max(muts) + 10])
    ax.set_ylabel(f'Count (n={len(muts)})', fontsize=LABEL_FONTSIZE)
    ax.set_xlabel('Branch length [# Mutations]', fontsize=LABEL_FONTSIZE)
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92, hspace=0.5)
    return fig


def plot_weights(data_in, out_file=None, bin_no=None):
    fig = plt.figure()
    gs = GridSpec(3, 1)

    y_labels = ['Topology FN', 'Avg. Mutation Prob', 'Topology FN+FP']
    if not bin_no:
        bin_no_best = int(np.sqrt(data_in.shape[0]))
        bin_no = np.linspace(data_in.min(), data_in.max(), bin_no_best)

    for i, data in enumerate(data_in.T):
        ax = fig.add_subplot(gs[i])
        ax.hist(data, bins=bin_no, rwidth=0.8)
        ax.set_ylabel(f'Frequency\n{y_labels[i]}', fontsize=LABEL_FONTSIZE)
        ax.set_ylim([None, data_in.shape[0] / 2])
    ax.set_xlabel('weights', fontsize=LABEL_FONTSIZE)

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.99, top=0.99, hspace=0.5)
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def generate_pval_plot(in_object, out_file=None, bin_no=None):
    if isinstance(in_object, str):
        df = pd.read_csv(in_object, sep='\t', engine='python', skipfooter=1)
        models = ['_'.join(i.split('_')[1:]) for i in df.columns[1::6]]
    elif isinstance(in_object, pd.DataFrame):
        df = in_object
        models = ['_'.join(i.split('_')[1:]) for i in df.columns[1::6]]
    else:
        raise IOError(f'Unknown input type: {type(in_object)}')

    fig = plt.figure(figsize=(16, 3 * len(models)))
    gs = GridSpec(3 * len(models), 1)
    for i, model in enumerate(models):
        ax = fig.add_subplot(gs[3*i:3*(i+1), 0])
        p_vals = df[f'p-value_{model}']
        if bin_no:
            bin_no_in = bin_no
        else:
            bin_no_in = max(15, int(np.sqrt(len(p_vals))))
        ax.hist(p_vals, bins=bin_no_in, range=(0, 1), alpha=0.5)
        ax.axhline(p_vals.size / bin_no_in, ls='--', c='red')
        ax.set_ylabel(f'counts {model}\n(n={p_vals.size})', fontsize=LABEL_FONTSIZE)
        ax.set_xlim([-0.01, 1.01])

    ax.set_xlabel('p-values', fontsize=LABEL_FONTSIZE)

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.99, top=0.99, hspace=0.5)
    
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def generate_gamma_plot(alphas, out_file):
    from scipy.stats import gamma
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


def generate_nbinom_plot(n, p, out_file):
    from scipy.stats import nbinom
    
    x = np.arange(nbinom.ppf(0.001, n, p), nbinom.ppf(0.999, n, p))
    y = nbinom.pmf(x, n, p)

    mu, var = nbinom.stats(n, p, moments='mv')
    half_stdev = np.sqrt(var/2)

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.plot(x, y, label=f'r={n:.2f},p={p:.2f}')

    ax.axvline(mu, ls='-', c='black')
    ax.axvline(mu - half_stdev, ls='--', c='black')
    ax.axvline(mu + half_stdev, ls='--', c='black')

    ax.set_ylabel(r'NB(k; r, p)', fontsize=LABEL_FONTSIZE)
    ax.set_ylim([0, y.max() * 1.05])
    ax.set_xlabel('k', fontsize=LABEL_FONTSIZE)
    ax.set_xlim([x.min() * 0.95, x.max() * 1.05])

    ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=TICK_FONTSIZE)
    ax.set_title('PDF', fontsize=LABEL_FONTSIZE * 1.5)

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_tree_matrix(X, out_file=None):
    fig, ax = plt.subplots(figsize=(12, 16))
    ax.matshow(X, cmap='Reds')
    ax.set_xlabel('Lambda', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Branch', fontsize=LABEL_FONTSIZE)
    ax.set_xticks([i for i in range(X.shape[1])])
    ax.set_xticklabels([str(i) for i in range(X.shape[1])])
    ax.set_yticks([i for i in range(X.shape[0])])
    ax.set_yticklabels([str(i) for i in range(X.shape[0])])

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()
    plt.close()


def plot_test_statistic(in_object, bin_no=None, out_file=None, in_dof=None):
    from scipy.stats.distributions import chi2

    if isinstance(in_object, str):
        df = pd.read_csv(in_object, sep='\t')
        models = ['_'.join(i.split('_')[1:]) for i in df.columns[1::6]]
    elif isinstance(in_object, pd.DataFrame):
        df = in_object
        models = ['_'.join(i.split('_')[1:]) for i in df.columns[1::6]]
    elif isinstance(in_object, np.ndarray):
        models = ['']
        df = pd.DataFrame(in_object, columns=['-2logLR_'])
        df['dof'] = in_dof
    else:
        raise IOError(f'Unknown input type: {type(in_object)}')

    fig = plt.figure(figsize=(16, 3 * len(models)))
    gs = GridSpec(3 * len(models), 1)

    max_x = 0
    for i, model in enumerate(models):
        ax = fig.add_subplot(gs[3*i:3*(i+1), 0])
        vals = df[f'-2logLR_{model}'].tolist()
        dof = df.iloc[-1][f'dof_{model}']
        if not bin_no:
            bin_no = max(15, int(np.sqrt(len(vals))))
        max_x = max(max_x, max(vals))
        ax.hist(vals, bins=bin_no, density=True , alpha=0.75)

        chi2_x = np.linspace(chi2.ppf(0.001, dof), chi2.ppf(0.999, dof), 1000)
        chi2_y = chi2.pdf(chi2_x, dof)
        ax.plot(chi2_x, chi2_y, 'r-', lw=5, alpha=0.8, label=r'$\chi^2$ (dof=' + f'{dof:.1f})')
            
        max_x = max(max_x, max(chi2_x))

        ax.set_ylabel(f'Density {model}', fontsize=LABEL_FONTSIZE)
        ax.set_xlim(0, int(max_x * 1.05))
        ax.legend(fontsize=TICK_FONTSIZE)
    ax.set_xlabel('Test statistic', fontsize=LABEL_FONTSIZE)

    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.99, top=0.99, hspace=0.5)

    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()

    # import pdb; pdb.set_trace()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input file')
    parser.add_argument('-f', '--format', type=str, default='pval',
        choices=['pval', 'mrbayes', 'gamma', 'nbinom', 'test'],
        help='What type of plot to generate')
    parser.add_argument('-a', '--alpha', type=float, default=[1], nargs='+',
        help='Alpha parameters of Gamma distribution with mean 1. Default = 1')
    parser.add_argument('-b', '--bins', type=int, default=0,
        help='Number of bins for histograms. Default = log(n)')
    parser.add_argument('-nbp', '--nbinom_pars', nargs=2, type=float,
        help='r und p parameters of negative binomial distribution.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
   
    if args.format == 'mrbayes':
        generate_mrbayes_plots(args.input,args.output)
    elif args.format == 'pval':
        generate_pval_plot(args.input, args.output, args.bins)
    elif args.format == 'test':
        plot_test_statistic(args.input, out_file=args.output, bin_no=args.bins)
    elif args.format == 'gamma':
        generate_gamma_plot(args.alpha, args.output)
    elif args.format == 'nbinom':
        generate_nbinom_plot(*args.nbinom_pars, args.output)
    else:
        raise IOError('Unknown format type: {}'.format(args.format))