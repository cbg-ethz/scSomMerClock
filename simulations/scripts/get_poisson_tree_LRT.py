#!/usr/bin/env python3

import os
import re
from io import StringIO

import numpy as np
# np.seterr(all='raise')
import pandas as pd
from scipy.stats.distributions import chi2
import scipy.special
from scipy.optimize import minimize, Bounds
from scipy.stats import poisson, nbinom
from scipy.spatial import distance

from utils import change_newick_tree_root
from plotting import _plot_pvals, plot_tree_matrix

from Bio import Phylo


import statsmodels.api as sm
from patsy import dmatrices

import warnings
from statsmodels.tools.sm_exceptions import DomainWarning
warnings.simplefilter('ignore', DomainWarning)
from statsmodels.genmod.families.links import identity

class pos_identity(identity):
    def __call__(self, p):
        return np.clip(p, 1e-3, 1e6)


LAMBDA_MIN = 1e-6
LAMBDA_MAX = np.inf


def simulate_poisson_tree(X, n=1000, pi_tresh=0):
    from tqdm import tqdm

    cols = [f"l{i}" for i in range(X.shape[1])]

    par_H0 = X.shape[1]
    X_H0_l = np.sum(X, axis=1) 

    bounds_H0 = np.full((par_H0, 2), (LAMBDA_MIN, np.inf))

    mu = np.clip(np.random.normal(100, 100, size=(n, par_H0)), LAMBDA_MIN, np.inf)
    # mu = np.clip(np.random.exponential(50, size=(n, par_H0)), LAMBDA_MIN, np.inf)

    # Zero inflated
    pi = np.random.random(size=(n, par_H0)) < pi_tresh
    mu = np.where(pi, 0, mu)

    X_mu = np.swapaxes((np.expand_dims(X, 1) * mu), 0, 1)
    Y = poisson.rvs(X_mu).sum(axis=2) + 1
    init_H0 = Y.mean(axis=1, keepdims=True) / X.sum(axis=0)

    mu_H0 = np.zeros((n, X.shape[0]))

    mu_H0_glm = np.zeros((n, X.shape[0]))

    for i in tqdm(range(n)):
        init_H0 = np.zeros(par_H0)
        for j, x_j in enumerate(X.T):
            # Lambda corresponding to exactly 1 count number
            rel = (x_j == 1) & (X[:,:j+1].sum(axis=1) == X_H0_l)
            rel_lgt = np.zeros(rel.sum())
            for k, rel_idx in enumerate(np.argwhere(rel).flatten()):
                rel_lgt[k] = max(LAMBDA_MIN, Y[i, rel_idx] \
                    - init_H0[np.argwhere(X[rel_idx][:j])].sum())
            init_H0[j] = rel_lgt.mean()

        opt_H0 = minimize(log_poisson, init_H0, args=(Y[i], X),
            bounds=bounds_H0)
        mu_H0[i] = np.dot(X, opt_H0.x)

        df = pd.DataFrame(X, columns=cols)
        df['Y'] = Y[i]
        expr = f'Y ~ {" + ".join(cols)} -1'
        y_glm, X_glm = dmatrices(expr, df, return_type='dataframe')
        glm = sm.GLM(y_glm, X_glm, family=sm.families.Poisson(pos_identity()))
        try:
            res_glm = glm.fit(start_params=init_H0,  maxiter=10000, disp=True)
        except:
            import pdb; pdb.set_trace()
            plot_tree_matrix(X)
        mu_H0_glm[i] = np.clip(np.dot(X, res_glm.params), LAMBDA_MIN, np.inf)

    ll_H0 = np.sum(poisson.logpmf(Y, mu_H0), axis=1)
    ll_H1 = np.sum(poisson.logpmf(Y, Y), axis=1)

    LR = -2 * (ll_H0 - ll_H1)
    p_vals = chi2.sf(LR, X.shape[0] - par_H0)

    ll_H0_glm = np.sum(poisson.logpmf(Y, mu_H0_glm), axis=1)
    LR_glm = -2 * (ll_H0_glm - ll_H1)
    p_vals_glm = chi2.sf(LR_glm, X.shape[0] - par_H0)

    show_pvals(p_vals)
    import pdb; pdb.set_trace()

    return p_vals


def simulate_nbinom_tree(X_H0, n=1000, pi_tresh=0):
    from tqdm import tqdm

    par_H1, par_H0 = X_H0.shape
    X_H1 = np.identity(par_H1)

    X_H0_l = np.sum(X_H0, axis=1) 

    bounds_H0 = Bounds((par_H0 + 1) * [LAMBDA_MIN], par_H0 * [np.inf] + [1 - LAMBDA_MIN])
    bounds_H0_short = Bounds((par_H0 + 0) * [LAMBDA_MIN], par_H0 * [np.inf] + [1 - LAMBDA_MIN])
    bounds_H1 = Bounds((par_H1 + 1) * [LAMBDA_MIN], par_H1 * [np.inf] + [1 - LAMBDA_MIN])
    bounds_H1_short = Bounds((par_H1 + 0) * [LAMBDA_MIN], par_H1 * [np.inf] + [1 - LAMBDA_MIN])

    mu = np.clip(np.random.normal(250, 100, size=(n, par_H0)), LAMBDA_MIN, np.inf)
    # mu = np.clip(np.random.exponential(50, size=(n, par_H0)), LAMBDA_MIN, np.inf)

    # Zero inflated
    pi = np.random.random(size=(n, par_H0)) < pi_tresh
    mu = np.where(pi, 0, mu)

    X_mu = np.clip(np.swapaxes((np.expand_dims(X_H0, 1) * mu), 0, 1), LAMBDA_MIN, np.inf)
    Y = np.zeros((n, par_H1))
    p = np.expand_dims(np.random.random(n), 1)
    p = np.clip(np.random.normal(0.5, 0.25), LAMBDA_MIN, 1 - LAMBDA_MIN)
    
    mu_H0 = np.zeros((n, par_H1))
    p_H0 = np.zeros((n, 1))
    mu_H1 = np.zeros((n, par_H1))
    p_H1 = np.zeros((n, 1))

    for i in tqdm(range(n)):
        Y[i] = nbinom.rvs(X_mu[i], p[i]).sum(axis=1)

        init_H0 = np.zeros(par_H0 + 1)
        for j, x_j in enumerate(X_H0.T):
            # Lambda corresponding to exactly 1 count number
            rel = (x_j == 1) & (X_H0[:,:j+1].sum(axis=1) == X_H0_l)
            rel_lgt = np.zeros(rel.sum())
            for k, rel_idx in enumerate(np.argwhere(rel).flatten()):
                rel_lgt[k] = max(0, Y[i, rel_idx] \
                    - init_H0[np.argwhere(X_H0[rel_idx][:j])].sum())
            init_H0[j] = rel_lgt.mean()
        init_H0[-1] = 0.5

        opt_H0 = minimize(log_nbinom, init_H0, args=(Y[i], X_H0), bounds=bounds_H0,
            options={'maxiter': 100000, 'maxfun': 1000000}, method='TNC')
        mu_H0[i] = np.dot(X_H0, opt_H0.x[:-1])
        p_H0[i] = opt_H0.x[-1]

        # opt_H0 = minimize(log_nbinom_short, init_H0[:-1], args=(Y[i], X_H0, p[i]),
        #     bounds=bounds_H0_short, options={'maxiter': 100000, 'maxfun': 1000000})
        # mu_H0[i] = np.dot(X_H0, opt_H0.x)
        # p_H0[i] = p[i]

        init_H1 = np.append(Y[i], 0.5)
        opt_H1 = minimize(log_nbinom, init_H1, args=(Y[i], X_H1), bounds=bounds_H1,
            options={'maxiter': 100000, 'maxfun': 1000000}, method='TNC')
        mu_H1[i] = np.dot(X_H1, opt_H1.x[:-1])
        p_H1[i] = opt_H1.x[-1]

        # opt_H1 = minimize(log_nbinom_short, Y[i], args=(Y[i], X_H1, p[i]),
        #     bounds=bounds_H1_short, options={'maxiter': 100000, 'maxfun': 1000000})
        # mu_H1[i] = np.dot(X_H1, opt_H1.x)
        # p_H1[i] = p[i]

    ll_H0 = np.sum(nbinom.logpmf(Y, mu_H0, p_H0), axis=1)
    ll_H1 = np.sum(nbinom.logpmf(Y, mu_H1, p_H1), axis=1)
    LR = -2 * (ll_H0 - ll_H1)
    p_vals = chi2.sf(LR, par_H1 - par_H0)

    show_pvals(p_vals)
    import pdb; pdb.set_trace()

    return p_vals


def get_mut_matrix(vcf_file, exclude_pat, include_pat):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude = []
    muts_raw = []
    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                    sample_no = len(sample_names)
                    samples = [0 for i in range(sample_no)]
                    if exclude_pat != '':
                        print('Exluded:')
                        for s_i, sample in enumerate(sample_names):
                            if re.fullmatch(exclude_pat, sample):
                                exclude.append(sample)
                                print('\t{}'.format(sample))
                        if len(exclude) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(exclude_pat))
                    if include_pat != '':
                        print('Include:')
                        for s_i, sample in enumerate(sample_names):
                            if not re.fullmatch(include_pat, sample):
                                exclude.append(sample)
                            else:
                                if not re.fullmatch(exclude_pat, sample):
                                    print('\t{}'.format(sample))
                        if len(exclude) == 0:
                            print('\nWARNING: no samples with pattern {} in vcf!\n' \
                                .format(include_pat))
                continue

            elif line.strip() == '':
                continue

            # VCF records
            line_cols = line.strip().split('\t')
            line_muts = np.zeros(sample_no)
            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')]
                # Missing in Monovar output format
                except ValueError:
                    line_muts[s_i] = np.nan
                    continue

                s_rec_ref = gt[0]
                s_rec_alt = gt[-1]
                if s_rec_ref == '.' or s_rec_alt == '.':
                    continue

                if s_rec_alt != '0' or s_rec_ref != '0':
                    line_muts[s_i] = 1
            muts_raw.append(line_muts)

    muts = pd.DataFrame(muts_raw, columns=sample_names)
    # Sanity check: remove sample without name
    exclude.append('')
    include = [i for i in sample_names if i not in exclude]
    
    return muts[include]


def write_tree(tree, out_file='example.newick'):
    tree.ladderize(reverse=False)
    Phylo.write(tree, out_file, 'newick')


def show_tree(tree):
    tree.ladderize(reverse=False)
    Phylo.draw(tree)


def get_tree_dict(tree_file, muts, paup_exe, min_dist=0):
    if 'cellphy_dir' in tree_file:
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True)
    elif 'trees_dir' in tree_file:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False)
    elif 'scite_dir' in tree_file:
        samples = [f'tumcell{i:0>4d}' for i in range(1, muts.shape[1], 1)] \
            + ['healthycell']
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            sample_names=samples, br_length=True)

    # With BioPython package
    tree = Phylo.read(StringIO(tree_str), 'newick')

    used_muts = []
    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_muts = muts[leaf_node.name]
        others_muts = muts.drop(leaf_node.name, axis=1).sum(axis=1)
        no_muts = ((leaf_muts == 1) & (others_muts == 0)).sum()
        used_muts.extend(muts[(leaf_muts == 1) & (others_muts == 0)].index.tolist())
        leaf_node.branch_length = no_muts
        
    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals():
        leaf_desc = [i.name for i in int_node.get_terminals()]
        leaf_no = len(leaf_desc)
        
        int_muts = muts[leaf_desc].sum(axis=1)
        others_muts = muts.drop(leaf_desc, axis=1).sum(axis=1)

        no_muts = ((int_muts == leaf_no) & (others_muts == 0)).sum()
        used_muts.extend(muts[(int_muts == leaf_no) & (others_muts == 0)].index.tolist())
        int_node.branch_length = no_muts
        int_node.name = '+'.join(leaf_desc)

    mut_diff = muts.shape[0] - tree.total_branch_length()
    if mut_diff > 0:
        print(f'Mutations that couldn\'t be mapped to the tree: {mut_diff}')
    # Collapse nodes with # mutations below threshold 
    write_tree(tree)
    while True:
        int_nodes = sorted([(tree.distance(i), i) \
                for i in tree.find_clades(terminal=False)],
            key=lambda x: x[0], reverse=True)
        # Iterate over internal nodes and check if criteria is met
        # if met, collapse node and redo whole procedure
        collapsed = False
        for _, int_node in int_nodes:
            br_length = [i.branch_length for i in int_node.clades]

            if sum(br_length) < min_dist or int_node.branch_length < min_dist:
                try:
                    print(f'\nCollapsing: {int_node.name}')
                    print(tree.total_branch_length())
                    parent = tree.collapse(int_node)
                    print(tree.total_branch_length())
                except ValueError:
                    pass
                else:
                    collapsed = True
                    break

        # No internal nodes fulfilled criteria
        if not collapsed:
            break
    write_tree(tree, 'tree_collapsed.newick')

    return tree


def get_ideal_tree(tree_file, mut_no=1000):
    # With BioPython package
    tree = Phylo.read(tree_file, 'newick')
    total_br = tree.total_branch_length()
    
    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_node.branch_length = leaf_node.branch_length / total_br * mut_no
        
    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals():
        try:
            int_node.branch_length = int_node.branch_length / total_br * mut_no
        except TypeError:
            int_node.name = '+'.join([i.name for i in tree.get_terminals()])
        else:
            int_node.name = '+'.join([i.name for i in int_node.get_terminals()])

    return tree


def get_model_data(tree):
    cell_no = tree.count_terminals()
    br_no = len([i for i in tree.find_clades()]) - 1
    br_indi_no = len([i for i in tree.find_clades(terminal=False)])

    X = np.zeros((br_no, br_indi_no), dtype=bool)
    Y = np.zeros(br_no, dtype=int)
    br_cells = {}
    
    br_dist = []
    for int_node in tree.find_clades(terminal=False):
        dist = tree.distance(int_node)
        no_desc = int_node.name.count('+')
        br_id = dist + (1 - no_desc / 1000)
        br_dist.append((br_id, int_node))
    sorted_br = sorted(br_dist, key=lambda x: x[0], reverse=True)

    # Iterate over internal nodes
    node_idx = 0
    for l_idx, (mut_no, int_node) in enumerate(sorted_br, 1):
        # Iterate over the two child nodes of each internal node
        for child in int_node.clades:
            # Set node properties
            Y[node_idx] = child.branch_length
            X[node_idx, 0:l_idx] = True
            # Subtract lambdas of child nodes
            for grand_child in child.clades:
                try:
                    X[node_idx] = X[node_idx] & ~X[br_cells[grand_child]]
                except:
                    import pdb; pdb.set_trace()
            # Sstore node index
            br_cells[child] = node_idx
            node_idx += 1
         
    # labels = [i[0].name for i in sorted(br_cells.items(), key=lambda x: x[1])]

    return Y, X.astype(int)


def log_poisson(x, k, T):
    return -np.nansum(poisson.logpmf(k, np.dot(T, x)))


def log_nbinom(x, k, T):
    return -np.nansum(nbinom.logpmf(k, np.dot(T, x[:-1]), x[-1]))


def log_nbinom_short(x, k, T, p):
    return -np.nansum(nbinom.logpmf(k, np.dot(T, x), p))


def opt_dist(x, k, T, dist_func=distance.minkowski, *args):
    return dist_func(np.dot(T, x), k, *args)


def show_pvals(p_vals):
    import matplotlib.pyplot as plt

    _plot_pvals(p_vals)
    plt.show()
    plt.close()


def _get_init(Y, X, p=None):
    X_l = np.sum(X, axis=1)

    init = np.zeros(X.shape[1])
    for j, x_j in enumerate(X.T):
        # Lambda corresponding to exactly 1 count number
        rel = (x_j == 1) & (X[:,:j+1].sum(axis=1) == X_l)
        rel_lgt = np.zeros(rel.sum())
        for k, rel_idx in enumerate(np.argwhere(rel).flatten()):
            rel_lgt[k] = max(0, Y[rel_idx] \
                - init[np.argwhere(X[rel_idx][:j])].sum())
        init[j] = rel_lgt.mean()
    if isinstance(p, float) and 0 < p < 1:
        return np.append(init, [p])
    else:
        return init


def get_LRT_poisson(Y, X, alpha=0.05):
    no_pars_H1, no_pars_H0 = X.shape

    init = _get_init(Y, X)

    bounds_H0 = np.full((no_pars_H0, 2), (LAMBDA_MIN, LAMBDA_MAX))
    bounds_H1 = np.full((X.shape[0], 2), (LAMBDA_MIN, LAMBDA_MAX))
    
    opt_H0 = minimize(log_poisson, init, args=(Y, X))#, bounds=bounds_H0)
    # opt_H0_dist = minimize(opt_dist, init, args=(Y, X), bounds=bounds_H0)
    # ll_H0_dist = np.sum(poisson.logpmf(Y, np.dot(X, opt_H0_dist.x)))

    import pdb; pdb.set_trace()
    ll_H0 = np.sum(poisson.logpmf(Y, np.dot(X, opt_H0.x)))
    ll_H1 = np.sum(poisson.logpmf(Y, Y))
    LR = -2 * (ll_H0 - ll_H1)
    # LR_dist = -2 * (ll_H0_dist - ll_H1)
    dof = no_pars_H1 - no_pars_H0
    p_val = chi2.sf(LR, dof)

    if p_val < alpha:
        hyp = 'H1'
    else:
        hyp = 'H0'

    return f'{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t{LR:0>5.2f}\t{p_val:.2E}\t{hyp}'


def get_glm_poisson_LRT(Y, X, alpha=0.05):
    init = _get_init(Y, X)
    
    glm = sm.GLM(pd.DataFrame(Y), pd.DataFrame(X),
        family=sm.families.Poisson(pos_identity()))
    res_glm = glm.fit(start_params=init,  maxiter=10000, disp=True)
    lambda_H0 = np.clip(np.dot(X, res_glm.params), LAMBDA_MIN, LAMBDA_MAX)

    # ll_H0 = res_glm.llf
    ll_H0 = np.sum(poisson.logpmf(Y, lambda_H0))
    ll_H1 = np.sum(poisson.logpmf(Y, Y))

    # LR = res_glm.deviance
    LR = -2 * (ll_H0 - ll_H1)

    dof = res_glm.df_resid
    p_val = chi2.sf(LR, dof)

    if p_val < alpha:
        hyp = 'H1'
    else:
        hyp = 'H0'

    return f'{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t{LR:0>5.2f}\t{p_val:.2E}\t{hyp}'


def get_LRT_nbinom(Y, X_H0, alpha=0.05):
    no_pars_H1, no_pars_H0 = X_H0.shape

    X_H0_l = np.sum(X_H0, axis=1)
    X_H1 = np.identity(no_pars_H1)

    bounds_H0 = np.append(np.full((no_pars_H0, 2), (LAMBDA_MIN, LAMBDA_MAX)),
        [[LAMBDA_MIN, 1 - LAMBDA_MIN]], axis=0)
    bounds_H1 = np.append(np.full((X_H0.shape[0], 2), (LAMBDA_MIN, LAMBDA_MAX)),
        [[LAMBDA_MIN, 1 - LAMBDA_MIN]], axis=0)

    init_H0 = _get_init(Y, X_H0, 0.5)
    init_H1 = np.append(Y, 0.5)

    opt_H0 = minimize(log_nbinom, init_H0, args=(Y, X_H0), bounds=bounds_H0,
        options={'maxiter': 100000, 'maxfun': 100000})

    opt_H1 = minimize(log_nbinom, init_H1, args=(Y, X_H1), bounds=bounds_H1,
        options={'maxiter': 100000, 'maxfun': 100000})

    
    ll_H0 = np.sum(nbinom.logpmf(Y, np.dot(X_H0, opt_H0.x[:-1]), opt_H0.x[-1]))
    ll_H1 = np.sum(nbinom.logpmf(Y, np.dot(X_H1, opt_H1.x[:-1]), opt_H1.x[-1]))
    LR = -2 * (ll_H0 - ll_H1)
    dof = no_pars_H1 - no_pars_H0
    p_val = chi2.sf(LR, dof)

    if p_val < alpha:
        hyp = 'H1'
    else:
        hyp = 'H0'

    return f'{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t{LR:0>5.2f}\t{p_val:.2E}\t{hyp}'


def test_poisson(vcf_file, tree_file, out_file, paup_exe, exclude='', include='',
            alpha=0.05):

    run = os.path.basename(vcf_file).split('.')[1]
    muts = get_mut_matrix(vcf_file, exclude, include)
    tree = get_tree_dict(tree_file, muts, paup_exe, 0)

    # tree = get_ideal_tree(tree_file)

    Y, X_H0 = get_model_data(tree)
    dof = X_H0.shape[0] - X_H0.shape[1]

    # Add 1 pseudocount
    Y += 1
    
    # simulate_nbinom_tree(X_H0, 1000, 0.0)
    # simulate_poisson_tree(X_H0, 1000, 0.25)

    models = [('poisson', get_LRT_poisson), ('nbinom', get_LRT_nbinom),
        ('poisson_glm', get_glm_poisson_LRT)]
    models = [('poisson', get_LRT_poisson), ('poisson_glm', get_glm_poisson_LRT)]
    # models = [('poisson_glm', get_glm_poisson_LRT)]

    cols = ['H0', 'H1', '-2logLR', 'p-value', 'hypothesis']
    header_str = 'run\tdof'
    model_str = f'{run}\t{dof}'
    for model_name, model_call in models:
        new_model_str = model_call(Y, X_H0, alpha)
        model_str += f'\t{new_model_str}'
        header_str += '\t' + '\t'.join([f'{i}_{model_name}' for i in cols])
        print(model_name, new_model_str)

    with open(out_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in vcf format')
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
        help='Significance threshold. Default = 0.05.')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''
        test_poisson(snakemake.input.vcf, snakemake.input.tree, 
            snakemake.output[0], snakemake.params.paup_exe,
            snakemake.params.exclude, snakemake.params.include)
    else:
        import argparse
        args = parse_args()
        test_poisson(args.vcf, args.tree, args.output, args.exe, args.exclude,
            args.include, args.alpha)