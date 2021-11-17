#!/usr/bin/env python3

import os
import re
import copy
import gzip
import tempfile
import subprocess
from io import StringIO

import numpy as np
np.seterr(all='raise')
np.set_printoptions(suppress=True)
import pandas as pd
from scipy.stats.distributions import chi2
from scipy.special import binom as binomCoeff
from scipy.optimize import minimize, Bounds, LinearConstraint
from scipy.stats import poisson, multinomial, betabinom, binom

from Bio import Phylo
from tqdm import tqdm
import nlopt

from utils import change_newick_tree_root


LAMBDA_MIN = 1e-6
log_LAMBDA_MIN = -50
MUT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def plot_distribution(d, b=500):
    import matplotlib.pyplot as plt

    d = d.flatten()
    d = d[~np.isnan(d)]

    fig, ax = plt.subplots(1, 1)
    ax.hist(d, histtype='stepfilled', bins=b)
    plt.show()
    plt.close()


def get_mut_df(vcf_file, exclude_pat, include_pat, filter=True):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude = []
    skip = 0
    data = [[], []]
    idx =[[], []]

    reads = []
    ref_gt = []
    with file_stream as f_in:
        for line in f_in:
            try:
                line = line.decode()
            except:
                pass
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
            # dims: 0 het, 1 hom, 2 true
            line_muts = np.zeros((3, sample_no))
            line_reads = np.zeros((sample_no, 4))
            diff_muts = []
            for s_i, s_rec in enumerate(line_cols[9:]):
                s_rec_elms = s_rec.split(':')
                try:
                    gt = s_rec_elms[0]
                    true_gt = s_rec_elms[-1]
                # Missing in Monovar output format
                except ValueError:
                    line_muts[:, s_i] = np.nan
                    continue

                # Get true gt
                if true_gt.count('0') == 1:
                    line_muts[2][s_i] = 1
                elif true_gt == '1|1' or true_gt == '2|2' or true_gt == '1|2' \
                        or true_gt == '2|1':
                    line_muts[2][s_i] = 2

                # Set filtered genotypes to nan
                if gt == '.|.':
                    line_muts[0:2, s_i] = np.nan
                elif gt == '0|0': # Wildtype
                    pass
                # Check if gt is mutation
                elif '0' in gt: # Heterozygous
                    line_muts[0][s_i] = int(gt[-1])
                    diff_muts.append(int(gt[-1]))
                else: # Homozygous
                    line_muts[1][s_i] = int(gt[-1])
                    diff_muts.append(int(gt[-1]))

                # Add reads
                try:
                    reads_i = np.array(s_rec_elms[2].split(','), dtype=int)
                except ValueError:
                    reads_i = np.full(4, np.nan)
                else:
                    if sum(reads_i) == 0:
                        reads_i = np.full(4, np.nan)

                line_reads[s_i] = reads_i

            # skip lines without any call, neither true or false
            if not any(np.nansum(line_muts, axis=1) > 0):
                skip += 1
                continue

            # Check if called mutations are all the same
            alts, alt_counts = np.unique(diff_muts, return_counts=True)

            if alts.size == 0:
                data[0].append(line_muts[0] + line_muts[1])
            # Two different single mutations: skip
            elif alts.size == 2 and all(alt_counts == 1):
                data[0].append(np.where(np.isnan(line_muts[0]), np.nan, 0))
                line_cols[6] = 'ISA_violation'
            else:
                ab_alt = alts[np.argmax(alt_counts)]
                het = np.where(line_muts[0] == ab_alt, 1, 0)
                hom = np.where(line_muts[1] == ab_alt, 2, 0)
                data[0].append(np.where(np.isnan(line_muts[0]), np.nan, hom+het))
                alt_idx = MUT[line_cols[4].split(',')[int(ab_alt) - 1]]

            data[1].append(line_muts[2])

            pos = int(line_cols[1])
            idx[1].append(pos)
            if 'PASS' in line_cols[6] or line_cols[6] == 's50':
                idx[0].append(pos)
                ref_gt.append((MUT[line_cols[3]], alt_idx))
                reads.append(line_reads)
            else:
                skip += 1

    muts = pd.DataFrame(data[0], index=idx[1], columns=sample_names)
    true_muts = pd.DataFrame(data[1], index=idx[1], columns=sample_names)

    # Sanity check: remove sample without name
    exclude.append('')
    include_id = [i for i in sample_names if i not in exclude]
    include_idx = [i for i,j in enumerate(sample_names) if j not in exclude]

    if filter:
        muts = muts.loc[idx[0], include_id]
        reads = np.array(reads)[:,include_idx]
        stats = get_stats(muts, true_muts.loc[idx[0], include_id], True, True)
    else:
        muts = muts[include_id]
        reads = np.array(reads)
        stats = get_stats(muts, true_muts, True, True)

    read_data = (np.array(idx[0]), np.array(ref_gt), reads, np.array(include_id))

    return muts, true_muts[include_id], read_data, stats


def get_stats(df, true_df, binary=True, verbose=False, per_cell=False):
    if binary:
        df = df.replace(2, 1)
        true_df = true_df.replace(2, 1)

    T = (df == true_df) & df.notnull()
    F = (df != true_df) & df.notnull()
    MS_df = df.isna()
    if per_cell:
        muts = T.sum()
        TP = (T & (df > 0)).sum() / muts
        TN = (T & (df == 0)).sum() / muts

        FP = (F & (df > true_df)).sum() / muts
        FN = (F & (df < true_df)).sum() / muts

        MS = MS_df.sum() / muts
        MS_N = (MS_df & (true_df == 0)).sum() / muts
        MS_T = (MS_df & (true_df > 0)).sum() / muts
    else:
        TP = (T & (df > 0)).sum().sum()
        TN = (T & (df == 0)).sum().sum()

        FP = (F & (df > true_df)).sum().sum()
        FN = (F & (df < true_df)).sum().sum()

        MS = MS_df.sum().sum()
        MS_N = (MS_df & (true_df == 0)).sum().sum()
        MS_T = (MS_df & (true_df > 0)).sum().sum()

    mut_wrong = (df[true_df.sum(axis=1) == 0].sum(axis=1) > 0).sum()
    mut_missing = (df[true_df.sum(axis=1) > 0].sum(axis=1) == 0).sum()
    if verbose and not per_cell:
        print(f'# Mutations: {df.shape[0]} ' \
            f'(wrong: {mut_wrong}; missing: {mut_missing})\n' \
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


def save_mut_mapping(tree, out_dir):
    tree = copy.deepcopy(tree)
    for i, node in enumerate(tree.find_clades()):
        # Rename to identify later
        if node.name.count('+') == 0:
            node.name = f'{i}.{node.name}'
        else:
            node.name = str(i)
        node.branch_length = node.mut_no_soft
        br_file = os.path.join(out_dir, f'muts_branch_{i:0>3}.tsv')
        with open(br_file, 'w') as f:
            f.write('\n'.join([f'1\t{i}' for i in node.muts_br]))
    write_tree(tree, out_file=os.path.join(out_dir, 'muts_branch_tree.newick'))


def show_tree(tree, dendro=False, br_length='mut_no_soft'):
    tree = copy.deepcopy(tree)
    tree.ladderize(reverse=False)

    if br_length == 'lambda':
        tree.root.lambd = 'root'
        br_labels = lambda c: c.lambd
    elif br_length == 'length':
        br_labels = lambda c: c.branch_length
    elif br_length == 'mut_no_hard':
        for i in tree.find_clades():
            i.branch_length = i.mut_no_hard
        try:
            tree.root.mut_no_true
            br_labels = lambda c: f'{c.mut_no_hard:.1f} ({c.mut_no_true})'
        except AttributeError:
            br_labels = lambda c: f'{c.mut_no_hard:.1f}'
    elif br_length == 'mut_no_soft':
        for i in tree.find_clades():
            i.branch_length = i.mut_no_soft
        try:
            tree.root.mut_no_true
            br_labels = lambda c: f'{c.mut_no_soft:.1f} ({c.mut_no_true})'
        except AttributeError:
            br_labels = lambda c: f'{c.mut_no_soft:.1f}'
    elif 'weights' in br_length:
        for i in tree.find_clades():
            i.branch_length = i.mut_no_soft

        try:
            tree.root.mut_no_true
        except AttributeError:
            true_muts_flag = False
        else:
            true_muts_flag = True

        if br_length == 'weights':
            br_labels = lambda c: f'{c.mut_no_soft:.0f} ({c.mut_no_true})\n' \
                + '/'.join([f'{i:.2f}' for i in c.weight_norm])
        elif br_length == 'weights_mut':
            if true_muts_flag:
                br_labels = lambda c: f'{c.mut_no_soft:.0f} ({c.mut_no_true})\n' \
                    f'{c.weight_norm[1][0]:.2f} +- {c.weight_norm[1][1]:.2f}'
            else:
                br_labels = lambda c: f'{c.mut_no_soft:.0f}\n' \
                    f'{c.weight_norm[1][0]:.2f} +- {c.weight_norm[1][1]:.2f}'
        elif br_length == 'mut_no_weights':
            if true_muts_flag:
                br_labels = lambda c: f'{c.mut_no_sp_count_distoft:.0f} ({c.mut_no_true})\n' \
                    + '/'.join([f'{i:.2f}' for i in c.weight_norm])
            else:
                br_labels = lambda c: f'{c.mut_no_soft:.0f}\n' \
                    + '/'.join([f'{i:.2f}' for i in c.weight_norm])
        elif br_length == 'weights_soft':
            if true_muts_flag:
                br_labels = lambda c: f'{c.mut_no_soft:.0f} ({c.mut_no_true})\n' \
                    f'{c.weight_norm[3]:.2f}'
            else:
                br_labels = lambda c: f'{c.mut_no_soft:.0f}\n' \
                    f'{c.weight["soft"]:.2f}'
    else:
        raise RuntimeError(f'Unknown branch length parameter: {br_length}')

    if dendro:
        max_depth = max(tree.depths().values())
        for i in tree.find_clades(terminal=False, order='postorder'):
            int_len = i.name.count('+')
            for child in i.clades:
                child.branch_length = int_len - child.name.count('+')

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


def get_xcoal_errors(FP, FN, c=1):
    FP = max(FP * c, LAMBDA_MIN)
    FN = max(FN * c, LAMBDA_MIN)
    TN_fin = np.log(1 - 2 * FP)
    TP_fin = np.log(1 - (FN * (2 * FP + 6 - 3 * FN) + 4 * FP) / 12)
    FP_fin = np.log(2 * FP)
    FN_fin = np.log((FN * (2 * FP + 6 - 3 * FN) + 4 * FP) / 12)
    return np.array([TP_fin, FP_fin, TN_fin, FN_fin]) # TP, FP, TN, FN


def get_scite_errors(FP, FN, c=1):
    FP = max(FP * c, LAMBDA_MIN) # alpha
    FN = max(FN * c, LAMBDA_MIN) # beta
    return np.log(np.array([1 - FN, FP, 1 - FP, FN])) # TP, FP, TN, FN


def get_tree(tree_file, paup_exe, samples=[], FN_fix=None, FP_fix=None,
            rates=False):
    fkt = 1
    FP = float(re.search('WGA0[\.\d,]*-0[\.\d]*-(0[\.\d]*)', tree_file).group(1)) \
        + 0.01
    # TODO <NB> hardcoded seq error. Take from config/directory structure
    FN = float(re.search('WGA(0[\.\d]*)[,\.\d]*?-', tree_file).group(1))
    errors = get_xcoal_errors(FP, FN, fkt)

    if 'cellphy' in tree_file:
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True,
            br_length=True)
        # tree_mapped, tree_approx = add_cellphy_mutation_map(tree_file, paup_exe, muts)
        # tree = tree_mapped
        log_file = tree_file.replace('.mutationMapTree', '.log')
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                log = f.read().strip()
        try:
            FP = float(re.search('SEQ_ERROR: (0.\d+(e-\d+)?)', log).group(1))
            FN = float(re.search('ADO_RATE: (0.\d+(e-\d+)?)', log).group(1))
        except AttributeError:
            pass
        else:
            errors = get_xcoal_errors(FP, FN, fkt)

    elif 'scite' in tree_file:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            sample_names=samples, br_length=True)
        # Get ML error rates from log file
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
        errors = get_scite_errors(FP, FN, fkt)
    else:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            br_length=True)

    if FP_fix and FN_fix:
        errors = get_scite_errors(FP_fix, FN_fix, fkt)

    tree = Phylo.read(StringIO(tree_str), 'newick')
    # Prune outgroup
    try:
        outg = [i for i in tree.find_clades() if i.name == 'healthycell'][0]
        tree.prune(outg)
    except IndexError:
        outg = None

    if rates:
        return tree, outg.name, (FN, FP)
    else:
        return tree, outg.name, errors


def get_tree_reads(tree_file, reads, paup_exe, FN_fix=None, FP_fix=None):
    tree, outg, errors = get_tree(tree_file, paup_exe, reads[3], FN_fix, FP_fix,
        True)
    M = map_mutations_CellCoal(tree, reads, errors, outg)

    MS = np.sum(np.sum(np.isnan(reads[2]), axis=2) > 0) \
        / (reads[2].shape[0] * reads[2].shape[1])
    FN = errors[0] + MS
    FP = errors[1]

    add_br_weigts(tree, np.log([1 - FN, FP, 1 - FP, FN]))

    return tree, FP, FN, M



def get_tree_gt(tree_file, muts, paup_exe, FN_fix=None, FP_fix=None):
    # Get rooted tree from newick file (outg removed) and errors
    samples = [f'tumcell{i:0>4d}' for i in range(1, muts.shape[1], 1)] \
        + ['healthycell']
    tree, outg, errors = get_tree(tree_file, paup_exe, samples, FN_fix, FP_fix)

    if outg in muts.columns:
        muts.drop(outg, axis=1, inplace=True)
    # Remove mutations that were only present in outgroup
    muts = muts[muts.sum(axis=1) > 0]

    M = map_mutations_Scite(tree, muts, errors)
    FP_ret = np.exp(errors)[1]
    FN_ret = np.exp(errors)[3]

    # For branch weighting: add missing rate to FN (errors: [TP, FP, TN, FN])
    MS = muts.isna().sum().sum() / muts.size
    errors[3] = np.log(np.exp(errors[3]) + MS) # FN
    errors[0] = np.log(1 - np.exp(errors[3])) # TN

    add_br_weigts(tree, errors)

    return tree, FP_ret, FN_ret, M


def get_rooted_tree(tree_str, paup_exe):
    temp_tree_file = tempfile.NamedTemporaryFile(delete=False)
    temp_tree_file.write(str.encode(tree_str))
    temp_tree_file.close()

    out_file = tempfile.NamedTemporaryFile(delete=False)
    paup_cmd = 'getTrees file={i};\n' \
        'DerootTrees;\n'\
        'outgroup healthycell;\n' \
        'RootTrees rootMethod=outgroup userBrLens=yes;\n' \
        'saveTrees format=Newick root=yes brLens=user taxaBlk=yes file={o} ;\n' \
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
    #     wrong_root_len = int(re.search(':(\d+)\):0;', tree_new).group(1))

    return Phylo.read(StringIO(tree_new), 'newick')


def add_cellphy_mutation_map(tree_file, paup_exe, muts):
    mut_file = tree_file.replace('mutationMapTree', 'mutationMapList')
    with open(mut_file, 'r') as f_muts:
        muts_raw = f_muts.read().strip().split('\n')

    with open(tree_file, 'r') as f_tree:
        tree_old = f_tree.read().strip()
    n = tree_old.count('cell')

    tree_old_approx = tree_old

    tree_init = Phylo.read(StringIO(tree_old), 'newick')
    tree_init.root_with_outgroup('healthycell')

    name_map = {}
    for node in tree_init.find_clades():
        terminals = '+'.join([i.name for i in node.get_terminals()])
        name_map[terminals] = node.name

    mut_map = {}
    for i in muts_raw:
        try:
            id, mut_no, muts = i.split('\t')
        except ValueError:
            id, mut_no = i.split('\t')
            muts = ''
        node_name = re.search(f'(\w*):\d+\.\d+\[{id}\]', tree_old).group(1)
        tree_old = re.sub(f'\d+\.\d+\[{id}\]', mut_no, tree_old)
        if muts:
            mut_map[node_name] = set([int(i[3:]) for i in muts.split(',')])
        else:
            mut_map[node_name] = set([])

    for l in re.findall('\d+\.\d+\[\d+\]', tree_old_approx):
        # TODO take value from config
        new_l = float(l.split('[')[0]) * 10000
        tree_old_approx = re.sub(re.escape(l), str(new_l), tree_old_approx)


    tree = get_rooted_tree(tree_old, paup_exe)
    for node in tree.find_clades():
        node.name = '+'.join([i.name for i in node.get_terminals()])
        node.mut_no_soft = node.branch_length

        try:
            node.muts_br = mut_map[name_map[node.name]]
        except KeyError:
            node.muts_br = set([])

        # Correct PAUP* rooting by replacing branch with 0 length branch
        # if node.name == 'healthycell':
        #     node.mut_no_soft = node.branch_length + wrong_root_len
        # elif node.name.count('+') == n - 2:
        #     assert node.branch_length == wrong_root_len, 'Cannot identify added root branch'
        #     node.mut_no_soft = 0
        #     node.branch_length = 0
        # else:
        #     node.mut_no_soft = node.branch_length

    tree_approx = get_rooted_tree(tree_old_approx, paup_exe)
    for node in tree_approx.find_clades():
        node.name = '+'.join([i.name for i in node.get_terminals()])
        node.mut_no_soft = node.branch_length

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


def _get_S(tree, cells, add_zeros=True):
    n = cells.size
    node_map = {}
    S = pd.DataFrame(data=np.zeros((2 * n - 1, n)), columns=cells)
    # Initialize node attributes on tree and get leaf nodes of each internal node
    for i, node in enumerate(tree.find_clades()):
        node_map[i] = node
        node.muts_br = set([])
        node.mut_no_soft = 0
        node.mut_no_hard = 0
        node.probs = []
        node.prob_weights = []
        leaf_nodes = [j.name for j in node.get_terminals()]
        S.loc[i, leaf_nodes] = 1
        node.name = '+'.join(leaf_nodes)

    S = S.values
    if add_zeros:
        S = np.vstack([S, np.zeros(n)])
    S_inv = 1 - S
    return S, S_inv, node_map


def map_mutations_Scite(tree, muts_in, errors):
    S, S_inv, node_map = _get_S(tree, muts_in.columns)

    idx_map = muts_in.index.values
    nans = np.isnan(muts_in).values
    # Convert ternary hetero/homozygous data to binary mutation data
    muts = np.where(nans, np.nan, muts_in.values.astype(bool).astype(float))
    muts_inv = 1 - muts

    # het = np.where(muts_in == 1, 1, 0)
    # het = np.where(nans, np.nan, het)
    # hom = np.where(muts_in == 2, 1, 0)
    # hom = np.where(nans, np.nan, hom)

    M = np.zeros((muts_in.shape[0], S.shape[0]))
    for i, mut in tqdm(enumerate(muts)):
        mut_data = np.stack([
            np.nansum(S * mut, axis=1), # TP
            np.nansum(S_inv * mut, axis=1), # FP
            np.nansum(S_inv * muts_inv[i], axis=1), # TN
            np.nansum(S * muts_inv[i], axis=1), # FN,
        ])
        probs = np.dot(errors, mut_data)
        probs_norm = _normalize_log_probs(probs)
        M[i] = probs_norm

        max_prob = probs_norm.max()
        best_nodes = np.argwhere(probs == np.max(probs)).flatten()

        for best_node in best_nodes:
            node_map[best_node].muts_br.add(idx_map[i])
            node_map[best_node].mut_no_hard += 1 / best_nodes.size
            node_map[best_node].probs.append(max_prob)
    # np.unique(np.argmax(M, axis=1), return_counts=True

    soft_assigned = M.sum(axis=0)
    for i, node in node_map.items():
        node.mut_no_soft = soft_assigned[i]

    return pd.DataFrame(M, index=muts_in.index)


def map_mutations_CellCoal(tree, read_data, errors, outg):
    # Remove outgroup mutations
    idx, gt, reads, cells = read_data
    if outg:
        cell_idx = np.argwhere(cells != outg).flatten()
        reads = reads[:,cell_idx,:]
        cells = cells[cell_idx]

    S, S_inv, node_map = _get_S(tree, cells)

    gamma = np.log(errors[0]/2)
    gamma_not = np.log(1 - errors[0])
    eps = np.log(errors[1] / 3)
    eps_not = np.log(1 - errors[1])

    eps_het = np.log(2/3 * errors[1])
    eps_het_not = np.log(1 - 2/3 * errors[1])
    eps_const = np.log(2)

    M = np.zeros((reads.shape[0], S.shape[0]))
    max_probs = np.zeros(reads.shape[0])

    for i in tqdm(range(idx.size)):
        reads_i = reads[i]
        gt_i = gt[i]
        depth = reads_i.sum(axis=1)

        # ref, !ref counts
        ref = reads_i[:,gt_i[0]]
        ref_not = depth - ref
        # alt, !alt counts
        alt = reads_i[:,gt_i[1]]
        alt_not = depth - alt
        # ref|alt, !(ref|alt) counts
        het = ref + alt
        het_not = depth - het
        # Get log prob p(b|0,0) for wt and p(b|1,1) homoyzgous genotypes
        # Clip to avoid that single cells outweights whole prob
        clip_min = -200
        p_ref = np.clip(ref * eps_not + ref_not * eps, clip_min, 0)
        p_alt = np.clip(alt * eps_not + alt_not * eps, clip_min, 0)
        # Get log pro for p(b|0,1) for het genotype
        noNAN = ~np.isnan(p_ref)
        p_count_dist = binom.logpmf(np.nan_to_num(alt), np.nan_to_num(depth), 0.5)

        p_noADO = het * eps_het_not + het_not * eps_het - depth * eps_const
        p_ADO = np.logaddexp(gamma + p_ref, gamma + p_alt, where=noNAN)

        p_het = np.logaddexp(gamma_not + p_noADO, p_ADO, where=noNAN)

        probs_full = np.clip(S_inv * p_ref + S * p_het + S * p_count_dist, clip_min, 0)
        probs_br = np.nansum(probs_full, axis=1)

        # plot_distribution(probs_full)
        probs_norm = _normalize_log_probs(probs_br)
        M[i] = probs_norm

        max_prob = probs_norm.max()
        best_nodes = np.argwhere(probs_br == np.max(probs_br)).flatten()
        max_probs[i] = max_prob

        # if idx[i] in [60, 67, 81, 110, 138, 263]: import pdb; pdb.set_trace()
        # if idx[i] == 3071: import pdb; pdb.set_trace()

        for best_node in sorted(best_nodes, reverse=True):
            # Remove muts where MLE contains full wt (outside of tree/all zero)
            try:
                node_map[best_node].muts_br.add(idx[i])
                node_map[best_node].mut_no_hard += 1 / best_nodes.size
                node_map[best_node].probs.append(max_prob)
            except KeyError:
                pass

    # Remove muts where MLE contains full wt (outside of tree/all zero)
    skip = np.where(M[:,-1] == M.max(axis=1), True, False)
    soft_assigned = M[~skip].sum(axis=0)

    for i, node in node_map.items():
        node.mut_no_soft = soft_assigned[i]

    return pd.DataFrame(M[~skip], index=idx[~skip])


def map_mutations_naive(tree, muts, min_dist):
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


def add_br_weigts(tree, errors):
    m = tree.count_terminals()
    n = 2 * m - 1

    # Get successor matrix
    nodes = [i for i in tree.find_clades()]

    S = np.zeros((n, m), dtype=int)
    Y = np.zeros(n)
    for i, node in enumerate(nodes):
        cells = [int(i[-4:]) - 1 for i in node.name.split('+')]
        if cells[0] == 25 and len(cells) == 1: print(i)
        S[i, cells] = 1
        Y[i] = node.mut_no_soft

    # Add zero line for wildtype probabilities
    S = np.vstack([S, np.zeros(m)])
    S_inv = 1 - S

    errors_rev = np.array([errors[0], errors[3], errors[2], errors[1]])

    w = np.zeros((n + 1, n + 1))

    for i, y in enumerate(S):
        data = np.stack([
            np.nansum(S * y, axis=1), # TP
            np.nansum(S_inv * y, axis=1), # FP
            np.nansum(S_inv * (1 - y), axis=1), # TN
            np.nansum(S * (1 - y), axis=1), # FN,
        ])
        probs = np.dot(errors, data)
        probs_norm = _normalize_log_probs(probs)
        probs_rev = np.dot(errors_rev, data)
        probs_norm_rev = _normalize_log_probs(probs_rev)

        # print(probs_norm.max())
        # if i == 28: import pdb; pdb.set_trace()

        if i < n:
            w[i] = probs_norm
        # prob for wiltype cell
        else:
            w[i] = probs_norm
            continue

        nodes[i].odds = probs[i] / max(LAMBDA_MIN, (1 - probs[i]))
        nodes[i].weight = {'top': probs_norm[i],
            'top_2': np.mean([probs_norm_rev[i], probs_norm[i]])}

        if sum(nodes[i].prob_weights) > 1:
            avg = np.average(nodes[i].probs, weights=nodes[i].prob_weights)
            std = np.sqrt(np.average((nodes[i].probs - avg) ** 2,
                weights=nodes[i].prob_weights))
            nodes[i].weight['mut'] = (avg, std)
        else:
            nodes[i].weight['mut'] = (probs_norm[i], np.nan)

    for i, j in enumerate(1 - np.abs(w.sum(axis=0) - 1)):
        if i == n:
            continue
        nodes[i].weight['soft'] = j


def add_true_muts(tree, df_true):
    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        leaf_node.true_muts_node = set(df_true.index[df_true[leaf_node.name] > 0])
        # leaf_node.false_muts = leaf_node.muts.difference(leaf_node.true_muts)

    for int_node in tree.get_nonterminals(order='postorder'):
        child0 = int_node.clades[0]
        child1 = int_node.clades[1]
        int_node.true_muts_node = child0.true_muts_node \
            .intersection(child1.true_muts_node)
        # child0.false_muts = int_node.true_muts_node.difference(int_node.true_muts)
        child0.true_muts_br = child0.true_muts_node \
            .difference(int_node.true_muts_node)
        child0.false_muts_br = child0.muts_br.difference(child0.true_muts_br)
        child0.missing_muts_br = child0.true_muts_br.difference(child0.muts_br)
        child0.mut_no_true = len(child0.true_muts_br)

        child1.true_muts_br = child1.true_muts_node \
            .difference(int_node.true_muts_node)
        child1.false_muts_br = child1.muts_br.difference(child1.true_muts_br)
        child1.missing_muts_br = child1.true_muts_br.difference(child1.muts_br)
        child1.mut_no_true = len(child1.true_muts_br)

    tree.root.true_muts_br = tree.root.true_muts_node
    tree.root.false_muts_br = tree.root.muts_br.difference(tree.root.true_muts_br)
    tree.root.missing_muts_br = tree.root.true_muts_br.difference(tree.root.muts_br)
    tree.root.mut_no_true = len(tree.root.true_muts_br)


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
                parents[-2].mut_no_true += node.mut_no_true
                tree.collapse(node)
            # Node closest to root
            elif len(parents) == 1 and len(tree.root.clades) == 1:
                tree.root.mut_no += node.mut_no
                tree.root.mut_no_true += node.mut_no_true
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


def get_model_data(tree, pseudo_mut=0, true_data=False, fkt=1):
    # Lambdas are assigned branchwise from the tips to the root, following the
    #   classical strict molecular clock

    internal_nodes = [i for i in tree.find_clades(terminal=False)]

    leaf_no = tree.count_terminals()
    br_no = 2 * leaf_no - 2

    # Domensions: 0 soft assignment; 1 hard assignment
    Y = np.zeros((2, br_no), dtype=float)
    init = np.zeros(br_no, dtype=float)
    weights = np.zeros((br_no, 4), dtype=float)
    leafs = np.zeros(br_no, dtype=float)
    constr = np.zeros((leaf_no - 1, br_no), dtype=int)

    # Get nodes in reverse depth first order
    br_dist = []
    for int_node in internal_nodes:
        dist = tree.distance(int_node)
        no_desc = int_node.name.count('+')
        br_id = dist + (1 - no_desc / 1000)
        br_dist.append((br_id, int_node))
    sorted_br = sorted(br_dist, key=lambda x: x[0], reverse=True)

    def get_traversal(start_node, tree):
        traverse = np.zeros(br_no, dtype=int)
        traverse[br_cells[start_node]] = 1
        node = start_node
        while True:
            # Exit case: terminal node
            if len(node.clades) == 0:
                break
            if node.clades[0].count_terminals() >= node.clades[1].count_terminals():
                traverse[br_cells[node.clades[0]]] = 1
                node = node.clades[0]
            else:
                traverse[br_cells[node.clades[1]]] = 1
                node = node.clades[1]
        return traverse

    # Iterate over internal nodes
    node_idx = 0
    br_cells = {}
    for l_idx, (_, int_node) in enumerate(sorted_br):
        is_terminal = [i.is_terminal() for i in int_node.clades]
        # Both terminal nodes
        if all(is_terminal):
            # Assign new lambda to terminals
            for i, terminal in enumerate(int_node.clades):
                br_cells[terminal] = node_idx
                if true_data:
                    Y[1, node_idx] = terminal.mut_no_true + pseudo_mut
                else:
                    Y[1, node_idx] = terminal.mut_no_hard + pseudo_mut
                Y[0, node_idx] = terminal.mut_no_soft + pseudo_mut
                weights[node_idx] = [j[0] if i == 'mut' else j \
                    for i,j in terminal.weight.items()]
                constr[l_idx, node_idx] = 1

                terminal.lambd = f'+{l_idx}'
                leafs[node_idx] = 1
                if i == 0:
                    init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)
                    node_idx += 1
                else:
                    constr[l_idx, node_idx] = -1

        # One internal, one terminal node
        elif any(is_terminal):
            # Assign new lambda to terminal
            terminal = int_node.clades[is_terminal.index(True)]

            br_cells[terminal] = node_idx
            if true_data:
                Y[1, node_idx] = terminal.mut_no_true + pseudo_mut
            else:
                Y[1, node_idx] = terminal.mut_no_hard + pseudo_mut
            Y[0, node_idx] = terminal.mut_no_soft + pseudo_mut
            weights[node_idx] = [j[0] if i == 'mut' else j \
                for i,j in terminal.weight.items()]
            constr[l_idx, node_idx] = 1
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)

            leafs[node_idx] = 1
            terminal.lambd = f'+{l_idx}'
            node_idx += 1

            # Assign lambda sum to internal
            internal = int_node.clades[is_terminal.index(False)]

            br_cells[internal] = node_idx
            if true_data:
                Y[1, node_idx] = internal.mut_no_true + pseudo_mut
            else:
                Y[1, node_idx] = internal.mut_no_hard + pseudo_mut
            Y[0, node_idx] = internal.mut_no_soft + pseudo_mut
            weights[node_idx] = [j[0] if i == 'mut' else j \
                for i,j in internal.weight.items()]
            constr[l_idx] = get_traversal(terminal, tree) - get_traversal(internal, tree)

            internal.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(constr[l_idx]) if j != 0]))

        # Both internal nodes
        else:
            shorter, longer = sorted(int_node.clades, key=lambda x: x.branch_length)
            # Assign new lambda to shorter branch
            br_cells[shorter] = node_idx
            if true_data:
                Y[1, node_idx] = shorter.mut_no_true + pseudo_mut
            else:
                Y[1, node_idx] = shorter.mut_no_hard + pseudo_mut
            Y[0, node_idx] = shorter.mut_no_soft + pseudo_mut
            weights[node_idx] = [j[0] if i == 'mut' else j \
                for i,j in shorter.weight.items()]
            constr[l_idx, node_idx] = 1
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)

            shorter.lambd = f'+{l_idx}'
            node_idx += 1
            # Assign lambda sum to longer
            br_cells[longer] = node_idx
            if true_data:
                Y[1, node_idx] = longer.mut_no_true + pseudo_mut
            else:
                Y[1, node_idx] = longer.mut_no_hard + pseudo_mut
            Y[0, node_idx] = longer.mut_no_soft + pseudo_mut
            weights[node_idx] = [j[0] if i == 'mut' else j \
                for i,j in longer.weight.items()]
            constr[l_idx] = get_traversal(shorter, tree) - get_traversal(longer, tree)

            longer.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(constr[l_idx]) if j != 0]))

        init_val = np.dot(constr[l_idx], init)

        # Make sure init is within bounds
        if init_val >= LAMBDA_MIN:
            init[node_idx] = init_val
        else:
            # If not, increase latested added lambda value
            init[node_idx - 1] += -init_val + (0.1)
            init[node_idx] = np.dot(constr[l_idx], init)
        node_idx += 1

    # Add data for last node
    # if true_data:
    #     Y[1, node_idx] = int_node.mut_no_true
    # else:
    #     Y[1, node_idx] = int_node.mut_no_hard + pseudo_mut
    # Y[0, node_idx] = int_node.mut_no_soft + pseudo_mut
    # weights[node_idx] = int_node.weight

    # shorter.lambd = f'+{l_idx + 1}'
    # init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)

    assert np.allclose(constr @ init, 0, atol=LAMBDA_MIN), \
        'Constraints not fulfilled for x_0'
    assert (init >= LAMBDA_MIN).all(), \
        f'Init value smaller than min. distance: {init.min()}'

    # Normalize weights
    weights = weights ** fkt
    weights_norm = weights.shape[0] * weights / weights.sum(axis=0)

    for node in tree.find_clades():
        i = np.abs(weights[:, 0] - node.weight['top']).argmin()
        node.weight_norm = weights_norm[i]

    # from plotting import plot_weights
    # plot_weights(weights_norm)


    if true_data:
        Y_out =  Y[1]
    else:
        Y_out = Y[0]

    return Y_out, constr, init, weights_norm[:,0], leafs


def get_LRT_poisson(Y, constr, init, weights=np.array([]), short=True,
            alg='trust-constr'):

    # weights = np.ones(Y.size)
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
                (Y * np.log(np.where(l > LAMBDA_MIN, l, np.nan)) - l) * weights) \
                / scale_fac
    else:
        Y = Y.round()
        init = init.round()
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

    # rel_constr = np.argwhere(constr[:,np.argwhere(weights).flatten()].sum(axis=1)) \
    #     .flatten()
    # constr = constr[rel_constr]

    const = [{'type': 'eq', 'fun': lambda x: np.matmul(constr, x)}]
    init = np.clip(init, LAMBDA_MIN, None)
    bounds = np.full((Y.size, 2), (LAMBDA_MIN, Y.sum()))

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
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    if np.sum(np.abs(opt.x - Y)) < 1e-6:
        print('WARNING: init values unchanged -> optimization likely failed.')

    if short:
        ll_H0 = np.nansum((Y * np.log(opt.x) - opt.x) * weights)
    else:
        ll_H0 = np.sum(poisson.logpmf(Y, np.clip(opt.x, LAMBDA_MIN, None))\
            * weights)

    dof = Y.size - constr.shape[0]
    LR = -2 * (ll_H0 - ll_H1)

    # U = fun_jac(opt.x, Y).reshape(Y.size, 1)
    # K = fun_hess(opt.x, Y)
    # np.fill_diagonal(K, np.clip(np.diag(K), LAMBDA_MIN, None))
    # K_inv =  np.zeros(K.shape)
    # np.fill_diagonal(K_inv, 1 / np.diag(K))
    # var = np.diag(np.sqrt(K_inv ** 2))

    # print(f' First Bartlett identity: {np.mean(fun_jac(opt.x, Y))}')
    # print('Second Bartlett identity: '
    #     f'{np.mean(fun_hess(opt.x, Y)) + np.var(fun_jac(opt.x, Y))}')

    on_bound = np.sum(opt.x <= LAMBDA_MIN ** 0.5)
    if on_bound > 0:
        dof_diff = np.arange(on_bound + 1)
        p_vals = np.clip(chi2.sf(LR, dof - dof_diff), 1e-100, 1)
        p_val_weights = np.where(
            np.isnan(p_vals), 0, binomCoeff(on_bound, dof_diff))
        p_val = np.average(np.nan_to_num(p_vals), weights=p_val_weights)
    else:
        p_val = chi2.sf(LR, dof)

    return ll_H0, ll_H1, LR, dof, on_bound, p_val


# def get_LRT_poisson_nlopt(Y, constr, init, weights=np.array([]), short=True):
#     if weights.size == 0:
#         weights = np.ones(Y.size)

#     scale_fac = 1

#     if short:
#         def f(x, grad):
#             if grad.size > 0:
#                 grad[:] = -(Y / x - 1) * weights / scale_fac
#             return -np.nansum((Y * np.log(x) - x) * weights) / scale_fac

#     else:
#         def f(x, grad):
#             if grad.size > 0:
#                 grad[:] = -(Y / x - 1) * weights / scale_fac
#             return -np.sum(poisson.logpmf(Y, x) * weights) / scale_fac

#     def c(result, x, grad):
#         if grad.size > 0:
#             grad[:] = -constr / scale_fac
#             result[:] = -constr @ x

#     opt = nlopt.opt(nlopt.LD_SLSQP, Y.size) # LD_AUGLAG, LD_SLSQP, LN_COBYLA, GN_ISRES
#     opt.set_min_objective(f)
#     opt.set_xtol_rel(LAMBDA_MIN)
#     opt.set_lower_bounds(np.full(Y.size, LAMBDA_MIN))
#     opt.set_upper_bounds(np.full(Y.size, Y.sum()))
#     opt.add_equality_mconstraint(c, np.full(constr.shape[0], LAMBDA_MIN))

#     xopt = opt.optimize(init)

#     if short:
#         ll_H1 = np.nansum((Y * np.log(np.where(Y > 0, Y, 1)) - Y) * weights)
#         ll_H0 = np.nansum((Y * np.log(xopt) - xopt) * weights)
#     else:
#         ll_H1 = np.sum(poisson.logpmf(Y, Y) * weights)
#         ll_H0 = np.sum(poisson.logpmf(Y, xopt) * weights)

#     dof = weights.sum() - constr.shape[0]
#     LR = -2 * (ll_H0 - ll_H1)
#     # LR_test = -2 * (np.nansum(Y * np.log(opt2.x) - opt2.x) - ll_H1)

#     on_bound = np.sum(xopt <= LAMBDA_MIN)
#     if on_bound > 0:
#         dof_diff = np.arange(on_bound + 1)
#         weights = binomCoeff(on_bound, dof_diff)
#         p_vals = np.clip(chi2.sf(LR, dof - dof_diff), 1e-100, 1)
#         p_val = np.average(p_vals, weights=weights)
#     else:
#         p_val = chi2.sf(LR, dof)

#     return ll_H0, ll_H1, LR, dof + on_bound, p_val


def run_poisson_tree_test(vcf_file, tree_file, out_file, paup_exe, weight_fkt,
        exclude='', include='', muts_per_branch=False):
    run = os.path.basename(vcf_file).split('.')[1]
    alpha = 0.05

    cols = ['H0', 'H1', '-2logLR', 'dof', 'p-value', 'hypothesis']
    header_str = 'run\tfiltered\ttrue_muts\tSNVs\tTP\tFP\tTN\tFN\tMS\tMS_T\t'
    header_str += '\t'.join([f'{col}_poissonTree_{weight_fkt}' for col in cols])
    model_str = ''

    filter_muts = [True]
    use_true_muts = [True, False]

    for filter_type in filter_muts:
        muts, true_muts, reads, stats = \
            get_mut_df(vcf_file, exclude, include, filter=filter_type)

        tree, FP, FN, M = get_tree_reads(tree_file, reads, paup_exe)
        add_true_muts(tree, true_muts)

        # tree1, FP, FN, M = get_tree_gt(tree_file, muts, paup_exe)
        # add_true_muts(tree1, true_muts)

        if muts_per_branch:
            out_dir = f'{tree_file}.muts_per_branch'
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            save_mut_mapping(tree, out_dir)
            print(f'Mutation to branch mapping written to: {out_dir}')
            exit()
        # tree = prune_leafs(tree)

        for muts_type in use_true_muts:
            if weight_fkt == 0:
                Y, constr, init, weights, leafs = get_model_data(tree,
                    pseudo_mut=0, true_data=muts_type, fkt=1)
                weights = np.ones(Y.size)
            else:
                Y, constr, init, weights, leafs = get_model_data(tree,
                    pseudo_mut=0, true_data=muts_type, fkt=weight_fkt)

            model_str += f'{run}\t{filter_type}\t{muts_type}\t{muts.shape[0]}\t' \
                f'{stats["TP"]}\t{stats["FP"]}\t{stats["TN"]}\t{stats["FN"]}\t' \
                f'{stats["MS"]}\t{stats["MS_T"]}'

            if np.any(Y != 0):
                ll_H0, ll_H1, LR, dof, on_bound, p_val = \
                    get_LRT_poisson(Y, constr, init, weights, short=True)

                # if p_val < 0.05: import pdb; pdb.set_trace()
                if np.isnan(ll_H0):
                    continue
                hyp = int(p_val < alpha)
                model_str += f'\t{ll_H0:0>5.2f}\t{ll_H1:0>5.2f}\t{LR:0>5.2f}\t' \
                    f'{dof+on_bound}\t{p_val:.2E}\tH{hyp}\n'

    with open(out_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in vcf format')
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    parser.add_argument('-w', '--weights', type=float, default=0,
        help='Faktor for branch weighting. If <1: no weighting. If >=1: weights ' \
            'to the power of the value are used.')
    parser.add_argument('-mpb', '--muts_per_branch', action='store_true',
        help='Write the mutation to branch mapping and exit.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''

        run_poisson_tree_test(
            vcf_file=snakemake.input.vcf,
            tree_file=snakemake.input.tree,
            out_file=snakemake.output[0],
            paup_exe=snakemake.params.paup_exe,
            weight_fkt=float(snakemake.wildcards.tree_weight),
            exclude=snakemake.params.exclude,
            include=snakemake.params.include
        )
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = os.path.join(os.path.dirname(args.vcf),
                'poisson_tree.LRT.tsv')
        run_poisson_tree_test(
            vcf_file=args.vcf,
            tree_file=args.tree,
            out_file=args.output,
            paup_exe=args.exe,
            weight_fkt=args.weights,
            exclude=args.exclude,
            include=args.include,
            muts_per_branch=args.muts_per_branch
        )
