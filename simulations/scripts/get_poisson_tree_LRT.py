#!/usr/bin/env python3

import os
import re
from io import StringIO

import numpy as np
# np.seterr(all='raise')
import pandas as pd
from scipy.stats.distributions import chi2
import scipy.special
from scipy.optimize import minimize
from scipy.stats import poisson


from utils import change_newick_tree_root

import statsmodels.api as sm
from newick import read, loads
from Bio import Phylo

EPSILON = 1E-12


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


def get_tree_dict(tree_file, muts, paup_exe):
    if 'cellphy_dir' in tree_file:
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True)
    elif 'trees_dir' in tree_file:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False)
    elif 'scite_dir' in tree_file:
        samples = [f'tumcell{i:0>4d}' for i in range(1, cells + 1, 1)]
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            sample_names=samples)

    # With BioPython package
    tree = Phylo.read(StringIO(tree_str), 'newick')

    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_muts = muts[leaf_node.name]
        others_muts = muts.drop(leaf_node.name, axis=1).sum(axis=1)
        no_muts = ((leaf_muts == 1) & (others_muts == 0)).sum()
        leaf_node.branch_length = no_muts
        
    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals():
        leaf_desc = [i.name for i in int_node.get_terminals()]
        leaf_no = len(leaf_desc)
        
        int_muts = muts[leaf_desc].sum(axis=1)
        others_muts = muts.drop(leaf_desc, axis=1).sum(axis=1)

        no_muts = ((int_muts == leaf_no) & (others_muts == 0)).sum()
        int_node.branch_length = no_muts
        int_node.name = '+'.join(leaf_desc)

    return tree


def get_model_data(tree, muts):
    total_muts = muts.shape[0]
    cell_no = tree.count_terminals()
    br_no = 2 * cell_no - 2
    br_indi_no = cell_no - 1

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
    for l_idx, (mut_no, int_node) in enumerate(sorted_br, 1):
        # Iterate over the two child nodes of each internal node
        for c_idx, child in enumerate(int_node.clades):
            # Get and store node index
            node_idx = (l_idx - 1) * 2 + c_idx
            br_cells[child] = node_idx
            # Set node properties
            Y[node_idx] = child.branch_length
            X[node_idx, 0:l_idx] = True
            # Subtract lambdas of child nodes
            for grand_child in child.clades:
                X[node_idx] = X[node_idx] & ~X[br_cells[grand_child]]
         
    labels = [i[0].name for i in sorted(br_cells.items(), key=lambda x: x[1])]
    return Y, X.astype(int), labels


def opt_log_poisson(k, l):
    # P(x=k) = l^k * e^(-l) / k!
    # log(P(x=k)) = k * log(l) - l [- log(k!)]
    
    y = x * np.log(l) - l
    return np.sum(y)


def log_poisson(x, k, T):
    return -np.nansum(poisson.logpmf(k, np.dot(T, x)))
    # l_tree = np.dot(T, x) + EPSILON
    # return -np.sum(k * np.log(l_tree) - l_tree) # - np.log(scipy.special.factorial(k))
    

def log_zero_infl_poisson(x, k, T):
    pi = x[-1]

    zero_flag = np.where(k == 0, 1, 0).astype(bool)
    l_tree = np.dot(T, x[:-1]) + EPSILON

    try:
        zero_ll = np.log(pi + (1 - pi) * np.exp(-l_tree[zero_flag]))
    except FloatingPointError:
        import pdb; pdb.set_trace()

    try:
        ll = np.log(1 - pi) + k[~zero_flag] * np.log(l_tree[~zero_flag]) - l_tree[~zero_flag]
    except FloatingPointError:
        import pdb; pdb.set_trace()

    return -np.sum(zero_ll) - np.sum(ll)


def test_poisson(vcf_file, tree_file, out_file, paup_exe, exclude='', include='',
            alpha=0.05):

    run = os.path.basename(vcf_file).split('.')[1]
    muts = get_mut_matrix(vcf_file, exclude, include)
    tree = get_tree_dict(tree_file, muts, paup_exe)
    Y, X, labels = get_model_data(tree, muts)

    lambda_min = 0.01
    lambda_max = Y.max()

    ## Dummy data

    # lambd = np.array([25, 50, 100, 150])
    # lambd = np.random.exponential(np.mean(Y) * 2, size=X.shape[1]).round()
    # Y = np.zeros(X.shape[0])
    # for i, cell in enumerate(X.astype(bool)):
    #     Y[i] = np.sum(np.random.poisson(lambd[cell.astype(bool)]))

    no_pars_unconst = Y.size

    log_poisson(Y, Y, np.identity(no_pars_unconst))
    bounds_unconst = [(lambda_min, lambda_max)] * no_pars_unconst

    Y = np.append(Y, np.mean(Y == 0))
    bounds_unconst.append((1e-6, 1 - 1e-6))

    opt_unconst = minimize(log_zero_infl_poisson, Y, args=(Y[:-1], np.identity(no_pars_unconst)),
            bounds=bounds_unconst)
    opt_unconst2 = minimize(log_poisson, Y[:-1], args=(Y[:-1], np.identity(no_pars_unconst)),
            bounds=bounds_unconst[:-1])

    no_pars_clock = X.shape[1]
    bounds_clock = [(lambda_min, lambda_max)] * no_pars_clock

    bounds_clock.append((1e-6, 1 - 1e-6))
    # Get starting position: mutation count on corresponding branches
    init_clock = np.zeros(no_pars_clock)
    for i in np.argwhere(X.sum(axis=1) == 1).flatten():
        l_i = np.argwhere(X[i] == 1)[0]
        init_clock[l_i] = Y[i]
    # If branch is not associated with single count, set it to mean 
    init_clock[np.argwhere(init_clock == 0)] = np.mean(Y)

    init_clock = np.append(init_clock, np.mean(Y == 0))
    
    opt_clock = minimize(log_zero_infl_poisson, init_clock, args=(Y[:-1], X), bounds=bounds_clock)
    opt_clock2 = minimize(log_poisson, init_clock[:-1], args=(Y[:-1], X), bounds=bounds_clock[:-1])

    LR = 2 * (opt_clock.fun - opt_unconst.fun)
    # LR2 =  np.sum((muts - mean_muts)**2) / mean_muts
    dof = no_pars_unconst - no_pars_clock
    p_val = chi2.sf(LR, dof)

    if p_val < alpha:
        hyp = 'H1'
    else:
        hyp = 'H0'

    import pdb; pdb.set_trace()

    with open(out_file, 'w') as f_out:
        f_out.write(
            'run\tH0\tH1\t-2logLR\tdof\tp-value\thypothesis\n' \
            f'{run}\t{opt_clock.fun:0>5.2f}\t{opt_unconst.fun:0>5.2f}\t' \
            f'{LR:0>5.2f}\t{dof}\t{p_val:.2E}\t{hyp}\n'
        )


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