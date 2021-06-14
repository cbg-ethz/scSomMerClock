#!/usr/bin/env python3

import os
import re
from io import StringIO

import numpy as np
np.seterr(all='raise')
np.set_printoptions(suppress=True)
import pandas as pd
from scipy.stats.distributions import chi2
import scipy.special
from scipy.optimize import minimize, Bounds, LinearConstraint
from scipy.stats import poisson, nbinom, multinomial, expon
from scipy.spatial import distance

from utils import change_newick_tree_root
from plotting import _plot_pvals, plot_tree_matrix, plot_test_statistic

from Bio import Phylo

# import statsmodels.api as sm
# from patsy import dmatrices
# import autograd.numpy as anp
# from autograd.scipy.special import gammaln
# from autograd import grad

# import nlopt

# import warnings
# from statsmodels.tools.sm_exceptions import DomainWarning
# warnings.simplefilter('ignore', DomainWarning)
# from statsmodels.genmod.families.links import identity


LAMBDA_MIN = 1e-6
LAMBDA_MAX = np.inf


def simulate_poisson_tree(X, n=1000, pi_tresh=0, glm=True):
    from tqdm import tqdm
    pi_tresh = 0.25
    par_H1, par_H0 = X.shape

    # mu = np.clip(np.random.normal(100, 100, size=(n, par_H0)), LAMBDA_MIN, LAMBDA_MAX)
    # mu = nbinom.rvs(np.arange(1, par_H0 * 10, 10) ** 1.05, 0.8, size=(n, par_H0))
    mu = np.random.exponential(np.arange(1, par_H0+1, 1) ** 1.5, size=(n, par_H0))

    # Zero inflated
    pi = np.random.random(size=(n, par_H0)) < pi_tresh
    mu = np.where(pi, 0, mu)

    X_mu = np.swapaxes((np.expand_dims(X, 1) * mu), 0, 1)
    Y = poisson.rvs(X_mu).sum(axis=2)
    # Y += 1

    dof = par_H1 - par_H0
    LR = np.zeros(n)
    LR_glm = np.zeros(n)
    for i in tqdm(range(n)):
        ll_H0, ll_H1 = get_LRT_poisson(Y[i], X, return_on_border=True)
        LR[i] = -2 * (ll_H0 - ll_H1)

        if glm:
            ll_H0_glm, ll_H1_glm = get_LRT_poisson_glm(Y[i], X)
            LR_glm[i] = -2 * (ll_H0_glm - ll_H1_glm)

    p_vals = chi2.sf(LR, dof)
    p_vals_glm = chi2.sf(LR_glm, dof)
    import pdb; pdb.set_trace()

    show_pvals(p_vals)
    plot_test_statistic(LR, in_dof=dof)
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
        mu_H0[i] = np.matmul(X_H0, opt_H0.x[:-1])
        p_H0[i] = opt_H0.x[-1]

        # opt_H0 = minimize(log_nbinom_short, init_H0[:-1], args=(Y[i], X_H0, p[i]),
        #     bounds=bounds_H0_short, options={'maxiter': 100000, 'maxfun': 1000000})
        # mu_H0[i] = np.matmul(X_H0, opt_H0.x)
        # p_H0[i] = p[i]

        init_H1 = np.append(Y[i], 0.5)
        opt_H1 = minimize(log_nbinom, init_H1, args=(Y[i], X_H1), bounds=bounds_H1,
            options={'maxiter': 100000, 'maxfun': 1000000}, method='TNC')
        mu_H1[i] = np.matmul(X_H1, opt_H1.x[:-1])
        p_H1[i] = opt_H1.x[-1]

        # opt_H1 = minimize(log_nbinom_short, Y[i], args=(Y[i], X_H1, p[i]),
        #     bounds=bounds_H1_short, options={'maxiter': 100000, 'maxfun': 1000000})
        # mu_H1[i] = np.matmul(X_H1, opt_H1.x)
        # p_H1[i] = p[i]

    ll_H0 = np.sum(nbinom.logpmf(Y, mu_H0, p_H0), axis=1)
    ll_H1 = np.sum(nbinom.logpmf(Y, mu_H1, p_H1), axis=1)
    LR = -2 * (ll_H0 - ll_H1)
    p_vals = chi2.sf(LR, par_H1 - par_H0)

    show_pvals(p_vals)
    plot_test_statistic(LR, in_dof=15)
    import pdb; pdb.set_trace()

    return p_vals


def get_mut_df(vcf_file, exclude_pat, include_pat):
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


def show_tree(tree, dendro=False, br_lambdas=True):
    tree.ladderize(reverse=False)

    if dendro:
        max_depth = max(tree.depths().values())
        for i in tree.find_clades(terminal=False, order='postorder'):
            int_len = i.name.count('+')
            for child in i.clades:
                child.branch_length = int_len - child.name.count('+')

    if br_lambdas:
        tree.root.lambd = 'root'
        br_labels = lambda c: c.lambd
    else:
        br_labels = lambda c: c.branch_length

    Phylo.draw(tree, branch_labels=br_labels,
        label_func=lambda c: c.name if c.name.count('+') == 0 else '',
        subplots_adjust=({'left': 0.01, 'bottom': 0.01, 'right': 0.99, 'top': 0.99})
    ) 


def get_tree_dict(tree_file, muts, paup_exe, min_dist=0):
    if 'cellphy_dir' in tree_file:
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True,
            br_length=True)
    elif 'trees_dir' in tree_file:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            br_length=True)
    elif 'scite_dir' in tree_file:
        samples = [f'tumcell{i:0>4d}' for i in range(1, muts.shape[1], 1)] \
            + ['healthycell']
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            sample_names=samples, br_length=True)

    # With BioPython package
    tree = Phylo.read(StringIO(tree_str), 'newick')

    # Get max age if branch lengths display age instead of mutation number
    max_age = muts.sum(axis=0).max()

    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_node.muts = set(muts.loc[muts[leaf_node.name] == 1].index)
        
    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals(order='postorder'):
        int_node.name = '+'.join([i.name for i in int_node.get_terminals()])

        child0 = int_node.clades[0]
        child1 = int_node.clades[1]

        int_node.muts = child0.muts.intersection(child1.muts)
        child0.age = max_age - len(child0.muts)
        child1.age = max_age - len(child1.muts)
        child0.mut_no = len(child0.muts.difference(int_node.muts))
        child1.mut_no = len(child1.muts.difference(int_node.muts))

    tree.root.mut_no = len(tree.root.muts)
    tree.root.age = max_age
    tree.root.lambd = 'root'

    # Collapse nodes with # mutations below threshold 
    if min_dist > 0:
        write_tree(tree)
        tree = collapse_branches(tree, min_dist)
        write_tree(tree, 'example.collapsed.newick')

    return tree


def collapse_branches(tree, min_dist):
    while True:
        int_nodes = sorted([(tree.distance(i), i) \
                for i in tree.find_clades(terminal=False)],
            key=lambda x: x[0], reverse=True)
        # Iterate over internal nodes and check if criteria is met
        # if met, collapse node and redo whole procedure
        collapsed = False
        for _, int_node in int_nodes:
            br_length = [i.mut_no for i in int_node.clades]

            if sum(br_length) < min_dist or int_node.mut_no < min_dist:
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
    return tree


def get_ideal_tree(tree_file, mut_no=1000, sample='poisson', min_dist=0):
    # With BioPython package
    tree = Phylo.read(tree_file, 'newick')
    total_br = tree.total_branch_length()

    if sample == 'multinomial':
        clades = [i for i in tree.find_clades() if i.branch_length]
        prob = [i.branch_length / total_br for i in clades]
        
        tree.root.branch_length = 0
        tree.root.name = '+'.join([i.name for i in tree.get_terminals()])

        for idx, cnt in enumerate(multinomial.rvs(mut_no, prob)):
            clades[idx].name = '+'.join([i.name for i in clades[idx].get_terminals()])
            clades[idx].branch_length = cnt
    elif sample == 'poisson':
        # Add number of mutations to terminal nodes
        for leaf_node in tree.get_terminals():
            leaf_node.branch_length = poisson \
                .rvs(leaf_node.branch_length / total_br * mut_no)
            
        # Add number of mutations to internal nodes
        for int_node in tree.get_nonterminals():
            try:
                int_node.branch_length = poisson \
                    .rvs(int_node.branch_length / total_br * mut_no)
            except TypeError:
                int_node.name = '+'.join([i.name for i in tree.get_terminals()])
            else:
                int_node.name = '+'.join([i.name for i in int_node.get_terminals()])
    elif sample == 'perfect':
        # Add number of mutations to terminal nodes
        for leaf_node in tree.get_terminals():
            leaf_node.branch_length = round(leaf_node.branch_length / total_br * mut_no)
            
        # Add number of mutations to internal nodes
        for int_node in tree.get_nonterminals():
            try:
                int_node.branch_length = round(int_node.branch_length / total_br * mut_no)
            except TypeError:
                int_node.name = '+'.join([i.name for i in tree.get_terminals()])
            else:
                int_node.name = '+'.join([i.name for i in int_node.get_terminals()])
    else:
        raise RuntimeError(f'Unknown sampling type: {sample}')

    # Collapse nodes with # mutations below threshold 
    if min_dist > 0:
        tree = collapse_branches(tree, min_dist)

    return tree


def get_sign(x):
    if x < 0:
        return '-'
    else:
        return '+'


def get_model0_data(tree, min_dist=0):
    # Lambdas are assigned branchwise from the tips to the root, following the 
    #   classical strict molecular clock

    internal_nodes = [i for i in tree.find_clades(terminal=False)]

    br_indi_no = tree.count_terminals() - 1 
    br_no = br_indi_no + len(internal_nodes)
    
    X = np.zeros((br_no, br_indi_no), dtype=int)
    Y = np.zeros(br_no, dtype=int)
    constr = np.zeros((br_indi_no, br_no), dtype=int)
    init = np.zeros(br_no)

    # Get nodes in reverse depth first order
    br_dist = []
    for int_node in internal_nodes:
        dist = tree.distance(int_node)
        no_desc = int_node.name.count('+')
        br_id = dist + (1 - no_desc / 1000)
        br_dist.append((br_id, int_node))
    sorted_br = sorted(br_dist, key=lambda x: x[0], reverse=True)

    def get_traversal(start_node, tree, X, br_cells):
        X_traverse = np.zeros(br_indi_no, dtype=int)
        traverse = start_node
        while True:
            # Exit case: terminal node
            try:
                child0 = X[br_cells[traverse.clades[0]]]
            except IndexError:
                break
            else:
                child1 = X[br_cells[traverse.clades[1]]]

            if sum(child1 == -1) > sum(child0 == -1):
                X_traverse += child0
                traverse = traverse.clades[0]
            else:
                X_traverse += child1
                traverse = traverse.clades[1]
        return X_traverse

    # Iterate over internal nodes
    node_idx = 0
    br_cells = {}
    lambdas = {}
    for l_idx, (mut_no, int_node) in enumerate(sorted_br):
        is_terminal = [i.is_terminal() for i in int_node.clades]
        # Both terminal nodes
        if all(is_terminal):
            # Assign new lambda to terminals
            for i, terminal in enumerate(int_node.clades):
                Y[node_idx] = max(min_dist, terminal.mut_no)
                X[node_idx, l_idx] = 1
                br_cells[terminal] = node_idx
                terminal.lambd = f'+{l_idx}'
                if i == 0:
                    lambdas[l_idx] = node_idx
                    init[node_idx] = Y[node_idx]
                    node_idx += 1
            
        # One internal, one terminal node
        elif any(is_terminal):
            # Assign new lambda to terminal
            terminal = int_node.clades[is_terminal.index(True)]
            Y[node_idx] = max(min_dist, terminal.mut_no)
            X[node_idx, l_idx] = 1
            br_cells[terminal] = node_idx
            terminal.lambd = f'+{l_idx}'
            lambdas[l_idx] = node_idx
            init[node_idx] = Y[node_idx]
            node_idx += 1
            # Assign lambda sum to internal
            internal = int_node.clades[is_terminal.index(False)]
            Y[node_idx] = max(min_dist, internal.mut_no)
            X[node_idx, l_idx] = 1
            X[node_idx] -= get_traversal(internal, tree, X, br_cells)
            br_cells[internal] = node_idx
            internal.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(X[node_idx]) if j != 0]))
            
        # Both internal nodes
        else:
            shorter, longer = sorted(int_node.clades, key=lambda x: x.branch_length)
            # Assign new lambda to shorter branch
            Y[node_idx] = max(min_dist, shorter.mut_no)
            X[node_idx, l_idx] = 1
            br_cells[shorter] = node_idx
            shorter.lambd = f'+{l_idx}'
            lambdas[l_idx] = node_idx
            init[node_idx] = Y[node_idx]
            node_idx += 1
            # Assign lambda sum to longer
            Y[node_idx] = max(min_dist, longer.mut_no)
            X[node_idx, l_idx] = 1
            # Add (shortest) lambda traversal over shorter branch
            X[node_idx] += get_traversal(shorter, tree, X, br_cells)
            # Subtract (shortest) lambda traversal over longer branch
            X[node_idx] -= get_traversal(longer, tree, X, br_cells)
            br_cells[longer] = node_idx
            longer.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(X[node_idx]) if j != 0]))

        constr_idx = node_idx // 2
        constr[constr_idx, node_idx] = 1
        for i, j in enumerate(X[node_idx]):
            if j:
                constr[node_idx // 2, i*2] = -j
        init_val = np.matmul(X[node_idx], init[::2])
        # Make sure init is within bounds
        if init_val >= 0:
            init[node_idx,] = init_val
        else:
            init[node_idx-1,] += -init_val + 1
            init[node_idx,] = np.matmul(X[node_idx], init[::2])
        node_idx += 1

    assert (constr @ init == 0).all(), 'Constraints not fulfilled for x_0'

    Y += 1
    # show_tree(tree, dendro=True)
    return Y, X, constr, init


def get_model1_data(tree):
    # Lambdas are assigned horizontally to each inter-node distance
    
    internal_nodes = [i for i in tree.find_clades(terminal=False)]

    br_indi_no = len(internal_nodes)
    br_no = br_indi_no + tree.count_terminals() - 1
    
    X = np.zeros((br_no, br_indi_no), dtype=bool)
    Y = np.zeros(br_no, dtype=int)
    
    # Get nodes in reverse depth first order
    br_dist = []
    for int_node in internal_nodes:
        dist = tree.distance(int_node)
        no_desc = int_node.name.count('+')
        br_id = dist + (1 - no_desc / 1000)
        br_dist.append((br_id, int_node))
    sorted_br = sorted(br_dist, key=lambda x: x[0], reverse=True)

    # Iterate over internal nodes
    node_idx = 0
    br_cells = {}
    for l_idx, (mut_no, int_node) in enumerate(sorted_br, 1):
        # Iterate over the two child nodes of each internal node
        for child in int_node.clades:
            # Set node properties
            Y[node_idx] = child.mut_no
            X[node_idx, 0:l_idx] = True
            # Subtract lambdas of child nodes
            for grand_child in child.clades:
                X[node_idx] = X[node_idx] & ~X[br_cells[grand_child]]
            # Store node index
            child.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(X[node_idx]) if j != 0]))
            br_cells[child] = node_idx
            node_idx += 1
    
    # show_tree(tree, dendro=False)
    # import pdb; pdb.set_Trace()
    return Y, X.astype(int)


def log_poisson(x, k, T):
    l = np.matmul(T, x)
    return -np.nansum(k * np.log(l) - l)
    # return -np.nansum(poisson.logpmf(k, l))


def log_multinomial(x, k, T):
    try:
        return -np.sum(k * np.log(np.matmul(T, x)))
    except:
        import pdb; pdb.set_trace()
    return -multinomial.logpmf(k, k.sum(), np.matmul(T, x))


def log_multinomial_H1(x, k):
    return -multinomial.logpmf(k, k.sum(), x)


def log_ZIP(x, k, T):
    pi = x[-1]
    mu = x[:-1]

    zero = np.log(pi + (1 - pi) * np.exp(-mu))
    pos = poisson.logpmf(k, np.matmul(T, mu)) + np.log(1 - pi)
    import pdb; pdb.set_trace()


def log_nbinom(x, k, T):
    return -np.nansum(nbinom.logpmf(k, np.matmul(T, x[:-1]), x[-1]))


def log_nbinom_short(x, k, T, p):
    return -np.nansum(nbinom.logpmf(k, np.matmul(T, x), p))


def opt_dist(x, k, T, dist_func=distance.minkowski, *args):
    return dist_func(np.matmul(T, x), k, *args)


def show_pvals(p_vals):
    import matplotlib.pyplot as plt

    _plot_pvals(p_vals)
    plt.show()
    plt.close()


def _get_init(Y, X, p=None, mu_min=LAMBDA_MIN):
    X_l = np.sum(X, axis=1)

    init = np.zeros(X.shape[1])
    for j, x_j in enumerate(X.T):
        # Lambda corresponding to exactly 1 count number
        rel = (x_j == 1) & (X[:,:j+1].sum(axis=1) == X_l)
        if any(rel):
            rel_lgt = np.zeros(rel.sum())
            for k, rel_idx in enumerate(np.argwhere(rel).flatten()):
                rel_lgt[k] = max(0, Y[rel_idx] \
                    - init[np.argwhere(X[rel_idx][:j])].sum())
            init[j] = rel_lgt.mean()
        else:
            init[j] = 0

    init = np.clip(init, mu_min, LAMBDA_MAX)
    if isinstance(p, float) and 0 < p < 1:
        init = np.append(init, [p])

    return init


def _get_init_clock(Y, X):
    init = np.full(X.shape[1], -1, dtype=float)
    for x in np.argwhere(np.sum(X == 1, axis=1) == 1).flatten():
        lambd = np.argwhere(X[x] == 1).flatten()[0]
        if init[lambd] < 0:
            init[lambd] = Y[x]
        else:
            init[lambd] = (Y[x] + init[lambd]) / 2

    # Fit constraints
    for i in np.argwhere(np.matmul(X, init) < 0):
        diff = (Y[i] - np.matmul(X[i], init))[0]
        lambd_change = np.argwhere(X[i][0] == 1).flatten()
        init[lambd_change] += diff / lambd_change.size
   
    return np.clip(init, LAMBDA_MIN, LAMBDA_MAX)


def get_LRT_exponential(Y, X, tree=False):
    def fun_H0(x, Y, X):
        # f(k|l) = l * e^(-l *t * k)
        # log_f(k|l) ~ log(l) - l * t * k
        l = x[-1]
        t_cum = anp.matmul(X, x[:-1])
        return -np.sum(l * t_cum * Y)

    init_H0 = _get_init(Y, X)
    init_H0 = np.append(init_H0, [1e-6])
    bounds_H0 = np.full((X.shape[1] + 1, 2), (LAMBDA_MIN, LAMBDA_MAX))
    # Get constraints prohibiting negative lambdas
    opt_H0 = minimize(fun_H0, init_H0, args=(Y, X), bounds=bounds_H0)

    l = opt_H0.x[-1]
    t = opt_H0.x[:-1]
    import pdb; pdb.set_trace()

    ll_H0 = np.sum(poisson.logpmf(Y, np.matmul(X, np.clip(opt_H0.x, LAMBDA_MIN, LAMBDA_MAX))))
    ll_H1 = np.sum(poisson.logpmf(Y, np.clip(Y, LAMBDA_MIN, LAMBDA_MAX)))

    import pdb; pdb.set_trace()

    dof_diff = np.sum(opt_H0.x == LAMBDA_MIN) - np.sum(Y == 0)
    return ll_H0, ll_H1, dof_diff



def get_LRT_multinomial(Y, X, tree=False):
    n = Y.sum()

    init = _get_init(Y, X)
    init2 = np.clip(init / np.matmul(X, init).sum(), LAMBDA_MIN, 1 - LAMBDA_MIN)

    # init = np.random.random(init.shape) + LAMBDA_MIN 
    # init /= np.matmul(X, init).sum()

    bounds = np.full((X.shape[1], 2), (LAMBDA_MIN, 1 - LAMBDA_MIN))
    def const_one(x, T):
        return 1 - np.matmul(T, x).sum()
    cons = ({'type': 'eq', 'fun': const_one, 'args': (X,)})

    from scipy.optimize import check_grad
    def fun(p, Y, X):
        return np.sum(Y * np.log(np.matmul(X, p)))

    def grad(p, Y, X):
        p_cum = np.matmul(X, p)
        grad = np.zeros(p.size)
        for i, p_i in enumerate(X.T):
            grad[i] = np.sum(np.where(p_i, Y, 0) / p_cum)
        return grad
    check_grad(fun, grad, init, Y, X)

    opt = minimize(log_multinomial, init, args=(Y, X), bounds=bounds,
        constraints=cons, options={'maxiter': 100000, 'disp': True})
    opt2 = minimize(log_multinomial, init2, args=(Y, X), bounds=bounds, 
        constraints=cons, options={'maxiter': 1000})
    # opt3 = minimize(log_multinomial, init, args=(Y, X), bounds=bounds, 
    #     constraints=cons, options={'maxiter': 1000, 'verbose': 3},
    #     hess=lambda x, y, z: np.zeros((X.shape[1], X.shape[1])), method='trust-constr')

    bounds_H1 = np.full((Y.size, 2), (LAMBDA_MIN, 1 - LAMBDA_MIN))
    cons_H1 = ({'type': 'eq', 'fun': lambda x: 1 - np.sum(x)})
    init_H1 = np.random.random(Y.size)
    init_H1 /= init_H1.sum()
    opt_H1 = minimize(log_multinomial_H1, init_H1, args=(Y,), bounds=bounds_H1, 
        constraints=cons_H1)

    ll_H0 = multinomial.logpmf(Y, n, np.matmul(X, opt.x))
    ll_H1 = multinomial.logpmf(Y, n, Y / Y.sum())

    return ll_H0, ll_H1


def get_LRT_nbinom(Y, X_H0, tree=False, alpha=0.05):
    no_pars_H1, no_pars_H0 = X_H0.shape

    X_H1 = np.identity(no_pars_H1)

    bounds_H0 = np.append(np.full((no_pars_H0, 2), (LAMBDA_MIN, LAMBDA_MAX)),
        [[LAMBDA_MIN, 1 - LAMBDA_MIN]], axis=0)
    bounds_H1 = np.append(np.full((X_H0.shape[0], 2), (LAMBDA_MIN, LAMBDA_MAX)),
        [[LAMBDA_MIN, 1 - LAMBDA_MIN]], axis=0)

    init_H0 = _get_init(Y, X_H0, p=0.5, mu_min=1)
    init_H1 = np.append(Y, 0.5)

    opt_H0 = minimize(log_nbinom, init_H0, args=(Y, X_H0), bounds=bounds_H0,
        options={'maxiter': 100000, 'maxfun': 100000})

    opt_H1 = minimize(log_nbinom, init_H1, args=(Y, X_H1), bounds=bounds_H1,
        options={'maxiter': 100000, 'maxfun': 100000})

    ll_H0 = np.sum(nbinom.logpmf(Y, np.matmul(X_H0, opt_H0.x[:-1]), opt_H0.x[-1]))
    ll_H1 = np.sum(nbinom.logpmf(Y, np.matmul(X_H1, opt_H1.x[:-1]), opt_H1.x[-1]))

    return ll_H0, ll_H1


def get_LRT_poisson_test(Y, X, constr, init):
    def fun_H0(l, Y):
        # weight = poisson.pmf(np.floor(l), l)
        weight = 1
        return -np.sum(poisson.logpmf(Y, l) * weight)

    bounds_H0 = np.full((X.shape[0], 2), (1, Y.sum()))    
    const = [{'type': 'eq', 'fun': lambda x: np.matmul(constr, x)}]
    # const = LinearConstraint(constr, np.full(X.shape[0], LAMBDA_MIN), np.full(X.shape[0], Y.max()),)

    opt_H0 = minimize(fun_H0, init, args=(Y,), constraints=const, bounds=bounds_H0,
        options={'disp': False, 'maxiter': 20000},
        method='SLSQP')

    l_H0 = opt_H0.x
    l_H1 = Y
    # w_H0 = poisson.pmf(np.floor(l_H0), l_H0)
    w_H0 = 1
    # w_H1 = poisson.pmf(Y, l_H1)
    w_H1 = 1

    ll_H0 = np.nansum(poisson.logpmf(Y, l_H0) * w_H0)
    ll_H1 = np.sum(poisson.logpmf(Y, l_H1) * w_H1)
    
    if not opt_H0.success:
        print('\nFAILED OPTIMIZATION\n')
        return np.nan, np.nan, np.nan

    dof = constr.shape[1] - constr.shape[0]
    LR = -2 * (ll_H0 - ll_H1)

    on_bound = np.sum(opt_H0.x <= 1)
    if on_bound > 0:
        coeffs = scipy.special.binom(on_bound, np.arange(on_bound + 1))
        pval_coeff = np.zeros(on_bound + 1)
        for i in range(on_bound+1):
            pval_coeff[i] = chi2.sf(LR, dof - i)
        p_val = np.sum(coeffs * pval_coeff) / np.sum(coeffs)
    else:
        p_val = chi2.sf(LR, dof)
    # dof = constr.shape[1] - constr.shape[0]
    return ll_H0, ll_H1, LR, dof, p_val


def get_LRT_poisson(Y, X, tree=False):

    def fun_H0(l, Y, X):
        # f(k|l) = l^k * e^(-l) / k! ~ l^k * e^(-l)
        # log_f(k|l) ~ k * log(l) - l
        lambdas_cum = np.matmul(X, l)
        weights = poisson.pmf(np.floor(lambdas_cum), lambdas_cum)
        
        # prob_short = -np.sum((Y * np.log(lambdas_cum) - lambdas_cum) * weights)
        # return prob_short

        prob_long = -np.sum(poisson.logpmf(Y, lambdas_cum) * weights)
        return prob_long

    init_H0 = _get_init_clock(Y, X)
    bounds_H0 = np.full((X.shape[1], 2), (LAMBDA_MIN, Y.sum()))
    # bounds_H0 = np.full((X.shape[1], 2), (-np.inf, LAMBDA_MAX))

    def const_fun(x):
        return np.matmul(X, x) - np.full(X.shape[0], LAMBDA_MIN)

    const = [{'type': 'ineq', 'fun': const_fun}]
    
    opt_H0 = minimize(fun_H0, init_H0, args=(Y, X),
        constraints=const, bounds=bounds_H0,
        options={'disp': True, 'maxiter': 20000})
    
    # # import pdb; pdb.set_trace()
    # # n = X.shape[1]
    # # opt = nlopt.opt(algorithm, n)
    # # opt.set_min_objective(fun_H0)
    # # opt.set_lower_bounds(LAMBDA_MIN)
    # # opt.add_inequality_constraint(const_fun, tol=0)

    # # Constraints for ‘trust-constr’ are defined as a single object
    # const = LinearConstraint(X, np.full(X.shape[0], LAMBDA_MIN),
    #     np.full(X.shape[0], LAMBDA_MAX))
    # opt_H0 = minimize(fun_H0, init_H0, args=(Y, X),
    #     hess= lambda x, v, z: np.zeros((X.shape[1], X.shape[1])),
    #     constraints=const, bounds=bounds_H0,
    #     options={'disp': True, 'maxiter': 10000, 'verbose': 1},
    #     method='trust-constr')

    # if any(opt_H0.x < 0):
    #     print('Negative lambda')
    #     print(np.argwhere(opt_H0.x < 0).flatten())
    #     show_tree(tree, br_lambdas=False)
    #     for i in tree.find_clades(order='postorder'):
    #         if i == tree.root:
    #             continue
    #         val_str = '\n'
    #         for j, lambdas in enumerate(i.lambd.split(' ')):
    #             val_str += get_sign(int(lambdas[1:]))
    #             val_str += f'{opt_H0.x[int(lambdas[1:])]:.0f} '
    #         i.lambd = i.lambd + val_str
    #     show_tree(tree, dendro=True)
    #     import pdb; pdb.set_trace()

    # init_H1 = Y
    # bounds_H1 = np.full((Y.size, 2), (LAMBDA_MIN, LAMBDA_MAX))
    # opt_H1 = minimize(fun_H1, init_H1, jac=grad(fun_H1), args=(Y,),
    #     bounds=bounds_H1)
    # ll_H1 = np.sum(poisson.logpmf(Y, opt_H1.x[:-1]))

    # l_H0_def = np.matmul(X, opt_H0_def.x)
    # w_H0_def = poisson.pmf(np.floor(l_H0_def), l_H0_def)
    # ll_H0_def = np.nansum(poisson.logpmf(Y, l_H0_def) * w_H0_def)

    l_H0 = np.matmul(X, opt_H0.x)
    l_H1 = Y
    w_H0 = poisson.pmf(np.floor(l_H0), l_H0)
    w_H1 = poisson.pmf(Y, l_H1)

    ll_H0 = np.nansum(poisson.logpmf(Y, l_H0) * w_H0)
    ll_H1 = np.sum(poisson.logpmf(Y, l_H1) * w_H1)
    
    if not opt_H0.success:
        print('FAILED OPTIMIZATION')
        import pdb; pdb.set_trace()

    if np.isnan(ll_H0): import pdb; pdb.set_trace()
    # no_border = sum(opt_H0.x <= LAMBDA_MIN)
    # LR = -2 * (ll_H0 - ll_H1)

    # chi2_dof = X.shape[0] - X.shape[1] - no_border

    # if any(np.matmul(X, opt_H0.x) < 0): import pdb; pdb.set_trace()
    dof = X.shape[0] - 1 * np.sum(Y == 0) - X.shape[1] - 1 * np.sum(opt_H0.x <= LAMBDA_MIN)
    # dof = X.shape[0] - X.shape[1]
    return ll_H0, ll_H1, dof




def test_poisson(vcf_file, tree_file, out_file, paup_exe, exclude='', include='',
            alpha=0.05):

    run = os.path.basename(vcf_file).split('.')[1]
    muts = get_mut_df(vcf_file, exclude, include)
    tree = get_tree_dict(tree_file, muts, paup_exe, 0)
    
    Y, X_H0, constr, init = get_model0_data(tree)
    # Y, X_H0 = get_model1_data(tree)

    # n = 1000
    # p_vals = np.zeros(n)
    # for i in range(n):
    #     tree_id = get_ideal_tree(tree_file, mut_no=muts.shape[0])
    #     Y, X_H0 = get_model1_data(tree_id)
    #     ll_H0, ll_H1, dof_diff = get_LRT_poisson(Y, X_H0, tree)
    #     p_vals[i] = chi2.sf(-2 * (ll_H0 - ll_H1) / 1, dof_poisson)
    # show_pvals(p_vals)
    # import pdb; pdb.set_trace()

    # Add 1 pseudocount
    # Y += 1
    # simulate_nbinom_tree(X_H0, 1000, 0.0)
    # simulate_poisson_tree(X_H0, 1000, 0.2, glm=False)

    models = [('poisson', get_LRT_poisson), ('multinomial', get_LRT_multinomial),
        ('nbinom', get_LRT_nbinom)]
    models = [('poisson', get_LRT_poisson), ('multinomial', get_LRT_multinomial)]
    models = [('poisson', get_LRT_poisson)]
    # models = [('exponential', get_LRT_exponential)]

    cols = ['H0', 'H1', '-2logLR', 'p-value', 'hypothesis']
    header_str = 'run\tdof'
    model_str = f'{run}'

    for model_name, model_call in models:
        # ll_H0, ll_H1, dof = model_call(Y, X_H0, tree)

        ll_H0, ll_H1, LR, dof, p_val = get_LRT_poisson_test(Y, X_H0, constr, init)
        if np.isnan(ll_H0):
            continue
        hyp = f'H{int(p_val < alpha)}'

        model_str += f'\t{dof}\t{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t{LR:0>5.2f}\t' \
            f'{p_val:.2E}\t{hyp}'
        header_str += '\t' + '\t'.join([f'{col}_{model_name}' for col in cols])

        # print(model_str)
        # x1, x2 = model_call(Y_id, X_H0_id)
        # chi2.sf(-2 * (x1 - x2), dof)
    # import pdb; pdb.set_trace()

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