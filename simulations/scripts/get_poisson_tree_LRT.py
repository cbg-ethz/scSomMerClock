#!/usr/bin/env python3

import os
import re
import copy
import tempfile
import subprocess
from io import StringIO

import numpy as np
np.seterr(all='raise')
np.set_printoptions(suppress=True)
import pandas as pd
from scipy.stats.distributions import chi2
import scipy.special
from scipy.optimize import minimize, Bounds, LinearConstraint
from scipy.stats import poisson, nbinom, multinomial

from Bio import Phylo
from tqdm import tqdm
import nlopt

from utils import change_newick_tree_root


LAMBDA_MIN = 1e-6
log_LAMBDA_MIN = np.log(LAMBDA_MIN)
LAMBDA_MAX = np.inf


def get_mut_df(vcf_file, exclude_pat, include_pat, filter=True):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude = []
    skip = 0
    data = [[], []]
    idx =[[], []]
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
            ref = line_cols[3]
            line_muts = np.zeros((2, sample_no))
            for s_i, s_rec in enumerate(line_cols[9:]):
                try:
                    gt = s_rec[:s_rec.index(':')].split('|')
                    true_gt = s_rec[-3:].split('|')
                # Missing in Monovar output format
                except ValueError:
                    line_muts[:, s_i] = np.nan
                    continue

                # Get true gt
                if true_gt[0] != '0' or true_gt[1] != '0':
                    line_muts[1][s_i] = 1
                # Set filtered genotypes to nan
                if gt[0] == '.' or gt[1] == '.':
                    line_muts[0][s_i] = np.nan
                # Check if gt is mutation
                elif gt[0] != '0' or gt[1] != '0':
                    line_muts[0][s_i] = 1

            for i in [0, 1]:
                data[i].append(line_muts[i])

            pos = int(line_cols[1])
            idx[1].append(pos)
            if 'PASS' in line_cols[6]:
                idx[0].append(pos)
            else:
                skip += 1

    muts = pd.DataFrame(data[0], index=idx[1], columns=sample_names)
    true_muts = pd.DataFrame(data[1], index=idx[1], columns=sample_names)

    # Sanity check: remove sample without name
    exclude.append('')
    include = [i for i in sample_names if i not in exclude]

    if filter:
        muts = muts.loc[idx[0], include]
        true_muts = true_muts.loc[idx[0], include]
    else:
        muts = muts[include]
        true_muts = true_muts[include]

    no_true = (true_muts.sum(axis=1) > 0).sum()
    no_false = muts.shape[0] - no_true
    stats = get_stats(muts, true_muts, True)

    return muts, true_muts, stats


def get_stats(df, true_df, verbose=False):
    TP = ((true_df == 1) & (df == 1)).sum().sum()
    FP = ((true_df == 0) & (df == 1)).sum().sum()
    TN = ((true_df == 0) & (df == 0)).sum().sum()
    FN = ((true_df == 1) & (df == 0)).sum().sum()
    MS = df.isna().sum().sum()
    MS_N = ((true_df == 0) & (df.isna())).sum().sum()
    MS_T = ((true_df == 1) & (df.isna())).sum().sum()

    if verbose:
        print(f'# Mutations: {df.shape[0]} ' \
            f'(wrong: {(true_df.sum(axis=1) == 0).sum()})\n' \
            f'\tTP: {TP: >7}\t({TP / df.size:.3f})\n' \
            f'\tFP: {FP: >7}\t({FP / df.size:.3f})\n' \
            f'\tTN: {TN: >7}\t({TN / df.size:.3f})\n' \
            f'\tFN: {FN: >7}\t({FN / df.size:.3f})\n' \
            f'\tMS: {MS: >7}\t({MS / df.size:.3f})\n' \
            f'\tMS_N: {MS_N: >5}\t({MS_N / df.size:.3f})\n' \
            f'\tMS_T: {MS_T: >5}\t({MS_T / df.size:.3f})')

    return {'TP': TP, 'FP': FP, 'TN': TN, 'FN': FN, 'MS': MS, 'MS_N': MS_N,
        'MS_T': MS_T}



def write_tree(tree, out_file='example.newick'):
    tree.ladderize(reverse=False)
    Phylo.write(tree, out_file, 'newick')


def show_tree(tree, dendro=False, br_length='mut_no'):
    tree = copy.deepcopy(tree)
    tree.ladderize(reverse=False)

    if dendro:
        max_depth = max(tree.depths().values())
        for i in tree.find_clades(terminal=False, order='postorder'):
            int_len = i.name.count('+')
            for child in i.clades:
                child.branch_length = int_len - child.name.count('+')

    if br_length == 'lambda':
        tree.root.lambd = 'root'
        br_labels = lambda c: c.lambd
    elif br_length == 'length':
        br_labels = lambda c: c.branch_length
    elif br_length == 'mut_no':
        for i in tree.find_clades():
            i.branch_length = i.mut_no

        # br_labels = lambda c: f'{c.mut_no:.0f}'
        try:
            tree.root.true_mut_no
            br_labels = lambda c: f' {c.mut_no:.1f} ({c.true_mut_no})'
        except AttributeError:
            br_labels = lambda c: f'{c.mut_no:.1f}'

    elif br_length == 'true_mut_no':
        for i in tree.find_clades():
            i.branch_length = i.true_mut_no
        br_labels = lambda c: f'{c.true_mut_no}'
    else:
        raise RuntimeError(f'Unknown branch length parameter: {br_length}')

    import matplotlib
    import matplotlib.pyplot as plt
    plt.rc('font', **{'size': 14})

    fig = plt.figure(figsize=(10, 20), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, branch_labels=br_labels,
        label_func=lambda c: c.name if c.name and c.name.count('+') == 0 else '',
        subplots_adjust=({'left': 0.01, 'bottom': 0.01, 'right': 0.99, 'top': 0.99}),
        axes=ax
    )


def get_tree(tree_file, muts, paup_exe, FN_fix=None, FP_fix=None):
    if 'cellphy' in tree_file:
        # _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True,
        #     br_length=True)
        # tree = Phylo.read(StringIO(tree_str), 'newick')
        tree_mapped, tree_approx = add_cellphy_mutation_map(tree_file, paup_exe, muts)
        tree = tree_approx
    else:
        if 'scite' in tree_file:
            samples = [f'tumcell{i:0>4d}' for i in range(1, muts.shape[1], 1)] \
                + ['healthycell']
            tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
                sample_names=samples, br_length=True)
            error_file = tree_file.replace('_ml0.newick', '.errors.csv')
            log_file =  tree_file.replace('_ml0.newick', '.log')
            if os.path.exists(error_file):
                with open(error_file, 'r') as f:
                    errors_raw = f.read().strip().split('\n')
                FP, FN = [float(i) for i in errors_raw[1].split(',')]
            elif os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_raw = f.read()
                FN = float(
                    re.search('best value for beta:\\\\t(\d.\d+(e-\d+)?)', log_raw) \
                        .group(1))
                FP = float(
                    re.search('best value for alpha:\\\\t(\d.\d+(e-\d+)?)', log_raw) \
                        .group(1))
            else:
                FP = 1e-3
                FN = 0.15
        else:
            tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
                br_length=True)
            FP = float(re.search('WGA0[\.\d]*-(0[\.\d]*)', tree_file).group(1))
            FN = float(re.search('WGA(0[\.\d]*)-', tree_file).group(1))

        if FN_fix:
            FN = FN_fix
        else:
            FN = max(FN, 0.01)

        if FP_fix:
            FP = FP_fix
        else:
            FP = max(FP, 1e-6)

        tree = Phylo.read(StringIO(tree_str), 'newick')
        map_mutations_to_tree(tree, muts.copy(), FP_rate=FP, FN_rate=FN)
        # map_mutations_to_tree_naive(tree, muts, min_dist)

    outg = [i for i in tree.find_clades() if i.name == 'healthycell'][0]
    tree.prune(outg)
    return tree, FP, FN


def add_cellphy_mutation_map(tree_file, paup_exe, muts):
    mut_file = tree_file.replace('mutationMapTree', 'mutationMapList')
    with open(mut_file, 'r') as f_muts:
        muts_raw = f_muts.read().strip().split('\n')

    with open(tree_file, 'r') as f_tree:
        tree_old = f_tree.read().strip()
    n = tree_old.count('cell')

    tree_old_approx = tree_old

    for i in muts_raw:
        try:
            id, mut_no, _ = i.split('\t')
        except ValueError:
            id, mut_no = i.split('\t')
        tree_old = re.sub(f'\d+\.\d+\[{id}\]', mut_no, tree_old)

    for l in re.findall('\d+\.\d+\[\d+\]', tree_old_approx):
        # TODO take value from config
        new_l = float(l.split('[')[0]) * 10000
        tree_old_approx = re.sub(re.escape(l), str(new_l), tree_old_approx)

    def get_rooted_tree(tree_str):
        temp_tree_file = tempfile.NamedTemporaryFile(delete=False)
        temp_tree_file.write(str.encode(tree_str))
        temp_tree_file.close()

        out_file = tempfile.NamedTemporaryFile(delete=False)
        paup_cmd = 'getTrees file={i};\n' \
            'DerootTrees;\n'\
            'outgroup healthycell;\n' \
            'RootTrees rootMethod=outgroup userBrLens=yes;\n' \
            'saveTrees format=Newick root=yes brLens=user file={o} ;\n' \
            'quit;'.format(i=temp_tree_file.name, o=out_file.name)

        paup_file = tempfile.NamedTemporaryFile(delete=False)
        paup_file.write(str.encode(paup_cmd))
        paup_file.close()

        shell_cmd = ' '.join([paup_exe, '-n', paup_file.name ])#, '>', '/dev/null'])
        paup = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = paup.communicate()
        paup.wait()

        assert stderr == b'', str(stdout) + '\n' + str(stderr)

        with open(out_file.name, 'r') as f_tree:
            tree_new = f_tree.read().strip()
        out_file.close()

        # try:
        #     wrong_root_len = int(re.search(':(\d+),healthycell', tree_new).group(1))
        # except AttributeError:
        #     import pdb; pdb.set_trace()
        #     wrong_root_len = int(re.search(':(\d+)\):0;', tree_new).group(1))

        return Phylo.read(StringIO(tree_new), 'newick')

    tree = get_rooted_tree(tree_old)
    for node in tree.find_clades():
        node.name = '+'.join([i.name for i in node.get_terminals()])
        node.mut_no = node.branch_length
        # Correct PAUP* rooting by replacing branch with 0 length branch
        # if node.name == 'healthycell':
        #     node.mut_no = node.branch_length + wrong_root_len
        # elif node.name.count('+') == n - 2:
        #     assert node.branch_length == wrong_root_len, 'Cannot identify added root branch'
        #     node.mut_no = 0
        #     node.branch_length = 0
        # else:
        #     node.mut_no = node.branch_length

    tree_approx = get_rooted_tree(tree_old_approx)
    for node in tree_approx.find_clades():
        node.name = '+'.join([i.name for i in node.get_terminals()])
        node.mut_no = round(node.branch_length)

    return tree, tree_approx


def _normalize_log_probs(probs):
    max_i = np.argmax(probs)
    try:
        exp_probs = np.exp(probs[np.arange(probs.size) != max_i] \
            - probs[max_i])
    except FloatingPointError:
        exp_probs = np.exp(
            np.clip(probs[np.arange(probs.size) != max_i] - probs[max_i],
                log_LAMBDA_MIN, 0)
        )
    probs_norm = probs - probs[max_i] - np.log1p(np.sum(exp_probs))
    return np.exp(np.clip(probs_norm, log_LAMBDA_MIN, 0))


def map_mutations_to_tree(tree, muts, FP_rate=1e-4, FN_rate=0.2):
    n = muts.shape[1]
    node_map = {}
    X = pd.DataFrame(data=np.zeros((2 * n - 1, n)), columns=muts.columns)
    # Initialize node attributes on tree and get leaf nodes of each internal node
    for i, node in enumerate(tree.find_clades()):
        node_map[i] = node
        node.muts = set([])
        node.mut_no = 0
        leaf_nodes = [j.name for j in node.get_terminals()]
        X.loc[i, leaf_nodes] = 1
        node.name = '+'.join(leaf_nodes)

    X = X.values
    idx_map = muts.index.values
    muts = muts.values
    X_inv = 1 - X
    muts_inv = 1 - muts

    FP_n = max(FP_rate, LAMBDA_MIN)
    FN_n = max(FN_rate, LAMBDA_MIN)
    TP = np.log(1 - FN_n)
    FP = np.log(FP_n)
    TN = np.log(1 - FP_n)
    FN = np.log(FN_n)

    errors = np.array([TP, FP, TN, FN])

    soft_assigned = np.zeros(X.shape[0])
    for i, mut in tqdm(enumerate(muts)):

        mut_data = np.stack([
            np.nansum(X * mut, axis=1), # TP
            np.nansum(X_inv * mut, axis=1), # FP
            np.nansum(X_inv * muts_inv[i], axis=1), # TN
            np.nansum(X * muts_inv[i], axis=1), # FN,
        ])
        probs = np.dot(errors, mut_data)
        best_nodes = np.argwhere(probs == np.max(probs)).flatten()

        for best_node in best_nodes:
            node_map[best_node].muts.add(idx_map[i])
            node_map[best_node].mut_no += 1 / best_nodes.size
        soft_assigned += _normalize_log_probs(probs)
    for i, node in node_map.items():
        node.mut_no_soft = soft_assigned[i]


def map_mutations_to_tree_naive(tree, muts, min_dist):
    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_node.muts = set(muts.loc[muts[leaf_node.name] == 1].index)

    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals(order='postorder'):
        int_node.name = '+'.join([i.name for i in int_node.get_terminals()])

        child0 = int_node.clades[0]
        child1 = int_node.clades[1]

        int_node.muts = child0.muts.intersection(child1.muts)
        child0.mut_no = len(child0.muts.difference(int_node.muts))
        child1.mut_no = len(child1.muts.difference(int_node.muts))

    tree.root.mut_no = len(tree.root.muts)
    tree.root.lambd = 'root'

    # Collapse nodes with # mutations below threshold
    if min_dist > 0:
        tree = collapse_branches(tree, min_dist)


def add_true_muts(tree, true_muts):
    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_node.true_muts = set(true_muts.loc[true_muts[leaf_node.name] == 1].index)

    for int_node in tree.get_nonterminals(order='postorder'):
        child0 = int_node.clades[0]
        child1 = int_node.clades[1]
        int_node.true_muts = child0.true_muts.intersection(child1.true_muts)
        child0.true_mut_no = len(child0.true_muts.difference(int_node.true_muts))
        child1.true_mut_no = len(child1.true_muts.difference(int_node.true_muts))

    tree.root.true_mut_no = len(tree.root.true_muts)


def prune_leafs(tree_in):
    tree = copy.deepcopy(tree_in)
    # Prune leaf nodes
    for leaf in tree.get_terminals():
        leaf_muts = leaf.mut_no
        parent = tree.collapse(leaf)
        continue
    # Merge unifurcating branches
    while True:
        for node in tree.find_clades(order='postorder'):
            parents = tree.get_path(node)
            if len(parents) > 1 and len(parents[-2].clades) == 1:
                parents[-2].mut_no += node.mut_no
                parents[-2].true_mut_no += node.true_mut_no
                tree.collapse(node)
            # Node closest to root
            elif len(parents) == 1 and len(tree.root.clades) == 1:
                tree.root.mut_no += node.mut_no
                tree.root.true_mut_no += node.true_mut_no
                tree.collapse(node)

        if tree.is_bifurcating():
            break

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


def get_model_data(tree, min_dist=0, true_data=False):
    # Lambdas are assigned branchwise from the tips to the root, following the
    #   classical strict molecular clock

    internal_nodes = [i for i in tree.find_clades(terminal=False)]

    br_indi_no = tree.count_terminals() - 1
    br_no = br_indi_no + len(internal_nodes)

    X = np.zeros((br_no, br_indi_no), dtype=int)
    Y = np.zeros(br_no, dtype=float)
    # distance from leaf nodes
    d = np.ones((br_no, 2), dtype=float)

    constr = np.zeros((br_indi_no, br_no), dtype=int)
    init = np.zeros(br_no)
    init_p = np.zeros(br_no)

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
    for l_idx, (_, int_node) in enumerate(sorted_br):
        is_terminal = [i.is_terminal() for i in int_node.clades]
        # Both terminal nodes
        if all(is_terminal):
            # Assign new lambda to terminals
            for i, terminal in enumerate(int_node.clades):
                if true_data:
                    Y[node_idx] = terminal.true_mut_no
                else:
                    Y[node_idx] = max(min_dist, round(terminal.mut_no))
                d[node_idx] = [1, np.nan]
                X[node_idx, l_idx] = 1
                br_cells[terminal] = node_idx
                terminal.lambd = f'+{l_idx}'
                if i == 0:
                    lambdas[l_idx] = node_idx
                    init[node_idx] = max(min_dist, Y[node_idx], 2 * LAMBDA_MIN)
                    node_idx += 1

        # One internal, one terminal node
        elif any(is_terminal):
            # Assign new lambda to terminal
            terminal = int_node.clades[is_terminal.index(True)]
            if true_data:
                Y[node_idx] = terminal.true_mut_no
            else:
                Y[node_idx] = max(min_dist, round(terminal.mut_no))
            d[node_idx] = [1, np.nan]
            X[node_idx, l_idx] = 1
            br_cells[terminal] = node_idx
            terminal.lambd = f'+{l_idx}'
            lambdas[l_idx] = node_idx
            init[node_idx] = max(min_dist, Y[node_idx], 2 * LAMBDA_MIN)
            node_idx += 1
            # Assign lambda sum to internal
            internal = int_node.clades[is_terminal.index(False)]
            if true_data:
                Y[node_idx] = internal.true_mut_no
            else:
                Y[node_idx] = max(min_dist, round(internal.mut_no))
            d[node_idx] = [internal.clades[0].name.count('+') + 1,
                internal.clades[1].name.count('+') + 1]
            X[node_idx, l_idx] = 1
            X[node_idx] -= get_traversal(internal, tree, X, br_cells)
            br_cells[internal] = node_idx
            internal.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(X[node_idx]) if j != 0]))

        # Both internal nodes
        else:
            shorter, longer = sorted(int_node.clades, key=lambda x: x.branch_length)
            # Assign new lambda to shorter branch
            if true_data:
                Y[node_idx] = shorter.true_mut_no
            else:
                Y[node_idx] = max(min_dist, round(shorter.mut_no))
            d[node_idx] = [shorter.clades[0].name.count('+') + 1,
                shorter.clades[1].name.count('+') + 1]
            X[node_idx, l_idx] = 1
            br_cells[shorter] = node_idx
            shorter.lambd = f'+{l_idx}'
            lambdas[l_idx] = node_idx
            init[node_idx] = max(min_dist, Y[node_idx], 2 * LAMBDA_MIN)
            node_idx += 1
            # Assign lambda sum to longer
            if true_data:
                Y[node_idx] = longer.true_mut_no
            else:
                Y[node_idx] = max(min_dist, round(longer.mut_no))
            d[node_idx] = [longer.clades[0].name.count('+') + 1,
                longer.clades[1].name.count('+') + 1]
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
        if init_val >= max(2 * LAMBDA_MIN, min_dist):
            init[node_idx,] = init_val
        else:
            init[node_idx-1,] += -init_val + (1 + min_dist)
            init[node_idx,] = np.matmul(X[node_idx], init[::2])
        node_idx += 1

    assert np.allclose(constr @ init, 0, atol=LAMBDA_MIN), \
        'Constraints not fulfilled for x_0'
    assert (init >= max(min_dist, 2 * LAMBDA_MIN)).all(), \
        f'Init value smaller than min. distance: {init.min()}'

    return Y, constr, init, d


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

    return Y, X.astype(int)


def log_poisson(x, k, T):
    l = np.matmul(T, x)
    return -np.nansum(k * np.log(l) - l)
    # return -np.nansum(poisson.logpmf(k, l))


def get_LRT_poisson(Y, constr, init, weights=np.array([]), short=True,
            alg='trust-constr'):
    if weights.size == 0:
        weights = np.ones(Y.size)

    if short:
        ll_H1 = np.sum((Y * np.log(np.where(Y > 0, Y, 1)) - Y) * weights)
        for i in [10000, 1000, 100, 10, 1]:
            scale_fac = (ll_H1 // i) * i
            if scale_fac != 0:
                break

        def fun_opt(l, Y):
            return -np.nansum(
                (Y * np.log(np.where(l>LAMBDA_MIN, l, np.nan)) - l) * weights) \
                / scale_fac
    else:
        ll_H1 = np.sum(poisson.logpmf(Y, Y) * weights)
        for i in [100, 10, 1]:
            scale_fac = (-ll_H1 // i) * i
            if scale_fac != 0:
                break
        def fun_opt(l, Y):
            l[:] = np.clip(l, LAMBDA_MIN, None)
            return -np.sum(poisson.logpmf(Y, l) * weights) / scale_fac

    def fun_jac(l, Y):
        return -(Y / np.clip(l, LAMBDA_MIN, None) - 1) * weights / scale_fac

    def fun_hess(l, Y):
        return np.identity(Y.size) * (Y / np.clip(l, LAMBDA_MIN, None)**2) * weights / scale_fac

    const = [{'type': 'eq', 'fun': lambda x: np.matmul(constr, x)}]
    init = np.clip(init, LAMBDA_MIN, None)
    bounds = np.full((constr.shape[1], 2), (LAMBDA_MIN, Y.sum()))

    if alg == 'trust-constr':
        opt = minimize(fun_opt, init, args=(Y,), constraints=const,
            bounds=bounds, method='trust-constr', jac=fun_jac, hess=fun_hess,
            options={'disp': False, 'maxiter': 10000, 'barrier_tol': 1e-5})
    elif alg == 'SLSQP':
        opt = minimize(fun_opt, init, args=(Y,), constraints=const,
            bounds=bounds, method='SLSQP', jac=fun_jac,
            options={'disp': False, 'maxiter': 10000,})
    else:
        opt = minimize(fun_opt, init, args=(Y,), constraints=const,
            bounds=bounds, method=alg, jac=fun_jac, hess=fun_hess)

    if not opt.success:
        print('\nFAILED POISSON OPTIMIZATION\n')
        opt = minimize(fun_opt, init, args=(Y,), constraints=const,
            bounds=bounds, method='SLSQP', jac=fun_jac,
            options={'disp': False, 'maxiter': 10000,})
        if not opt.success:
            print('\nFAILED POISSON OPTIMIZATION, TWICE\n')
            return np.nan, np.nan, np.nan, np.nan, np.nan,

    if short:
        ll_H0 = np.nansum((Y * np.log(opt.x) - opt.x) * weights)
    else:
        ll_H0 = np.sum(poisson.logpmf(Y, np.clip(opt.x, LAMBDA_MIN, None))\
            * weights)

    dof = Y.size - constr.shape[0]
    LR = -2 * (ll_H0 - ll_H1)
    # LR_test = -2 * (np.nansum(Y * np.log(opt2.x) - opt2.x) - ll_H1)

    on_bound = np.sum(opt.x <= 10 * LAMBDA_MIN)
    if on_bound > 0:
        dof_diff = np.arange(on_bound + 1)
        p_vals = np.clip(chi2.sf(LR, dof - dof_diff), 1e-100, 1)
        p_val_weights = np.where(np.isnan(p_vals), 0,
            scipy.special.binom(on_bound, dof_diff))
        p_val = np.average(np.nan_to_num(p_vals), weights=p_val_weights)
    else:
        p_val = chi2.sf(LR, dof)

    return ll_H0, ll_H1, LR, dof + on_bound, p_val


# def get_LRT_poisson_nlopt(Y, X, constr, init, short=True):
#     scale_fac = 100

#     def f(x, grad):
#         if grad.size > 0:
#             grad[:] = -(Y / x - 1) / scale_fac
#         return -np.sum(poisson.logpmf(Y, x)) / scale_fac

#     def c(x, grad, c):
#         if grad.size > 0:
#             grad[:] = -c / scale_fac
#         return -np.sum(c * x) / scale_fac

#     def c1(result, x, grad):
#         if grad.size > 0:
#            grad[:] = -constr / scale_fac
#            result[:] = -constr @ x  / scale_fac


#     opt = nlopt.opt(nlopt.LD_SLSQP, Y.size) # LD_AUGLAG, LD_SLSQP, LN_COBYLA, GN_ISRES
#     opt.set_min_objective(f)
#     opt.set_xtol_rel(LAMBDA_MIN)
#     opt.set_lower_bounds(np.full(Y.size, LAMBDA_MIN))
#     opt.set_upper_bounds(np.full(Y.size, Y.sum()))
#     opt.add_equality_mconstraint(c1, np.full(constr.shape[0], 2 * LAMBDA_MIN))

#     xopt = opt.optimize(init)

#     ll_H1 = np.sum(poisson.logpmf(Y, Y))
#     ll_H0 = np.sum(poisson.logpmf(Y, xopt))

#     dof = Y.size - constr.shape[0]
#     LR = -2 * (ll_H0 - ll_H1)
#     # LR_test = -2 * (np.nansum(Y * np.log(opt2.x) - opt2.x) - ll_H1)

#     on_bound = np.sum(xopt <= 2 * LAMBDA_MIN)
#     if on_bound > 0:
#         dof_diff = np.arange(on_bound + 1)
#         weights = scipy.special.binom(on_bound, dof_diff)
#         p_vals = np.clip(chi2.sf(LR, dof - dof_diff), 1e-100, 1)
#         p_val = np.average(p_vals, weights=weights)
#     else:
#         p_val = chi2.sf(LR, dof)

#     return ll_H0, ll_H1, LR, dof + on_bound, p_val


def test_data(vcf_file, tree_file, out_file, paup_exe, exclude='', include='',
            weight=1):

    run = os.path.basename(vcf_file).split('.')[1]
    alpha = 0.05

    cols = ['H0', 'H1', '-2logLR', 'dof', 'p-value', 'hypothesis']
    header_str = 'run\tfiltered\ttrue_muts\tSNVs\tTP\tFP\tTN\tFN\tMS\tMS_T\t'
    header_str += '\t'.join([f'{col}_poissonTree' for col in cols])
    model_str = ''

    filter_muts = [False, True]
    use_true_muts = [True, False]

    for filter_type in filter_muts:
        muts, true_muts, stats = get_mut_df(vcf_file, exclude, include,
            filter=filter_type)
        tree, FP, FN = get_tree(tree_file, muts, paup_exe,)
            # FP_fix=stats['FP'] / muts.size,
            # FN_fix=(stats['FN'] + stats['MS_T']) / muts.size)
        add_true_muts(tree, true_muts)

        # tree = prune_leafs(tree)

        for muts_type in use_true_muts:
            Y, constr, init, leaf_dist = get_model_data(tree, min_dist=0,
                true_data=muts_type)
            weight_loss = FN
            weights = 1 - np.nansum(FN ** leaf_dist, axis=1) \
                - np.nansum(FP ** leaf_dist, axis=1)

            model_str += f'{run}\t{filter_type}\t{muts_type}\t{muts.shape[0]}\t' \
                f'{stats["TP"]}\t{stats["FP"]}\t{stats["TN"]}\t{stats["FN"]}\t' \
                f'{stats["MS"]}\t{stats["MS_T"]}'

            if np.any(Y != 0):
                ll_H0, ll_H1, LR, dof, p_val = get_LRT_poisson(Y, constr, init, weights)
                if np.isnan(ll_H0):
                    continue
                hyp = f'H{int(p_val < alpha)}'

                model_str += f'\t{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t{LR:0>5.2f}\t{dof}\t' \
                    f'{p_val:.2E}\t{hyp}\n'

    # if p_val < alpha: import pdb; pdb.set_trace()
    with open(out_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in vcf format')
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
    parser.add_argument('-w', '--weights', type=float, default=1,
        help='Weight for leaf branches (<1 = down-weighting). Default = 1')
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
        test_data(snakemake.input.vcf, snakemake.input.tree, 
            snakemake.output[0], snakemake.params.paup_exe,
            snakemake.params.exclude, snakemake.params.include,
            snakemake.params.weights)
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = os.path.join(os.path.dirname(args.vcf),
                'poisson_tree.LRT.tsv')
        test_data(args.vcf, args.tree, args.output, args.exe, args.exclude,
            args.include, args.weights)