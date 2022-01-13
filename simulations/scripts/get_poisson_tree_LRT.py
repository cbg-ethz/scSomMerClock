#!/usr/bin/env python3

import os
import re
import gzip
from io import StringIO

import numpy as np
np.seterr(all='raise')
np.set_printoptions(suppress=True)
import pandas as pd
from scipy.stats.distributions import chi2
from scipy.special import binom as binomCoeff
from scipy.optimize import minimize

from Bio import Phylo
from tqdm import tqdm

from utils import change_newick_tree_root


LAMBDA_MIN = 1e-6
log_LAMBDA_MIN = -50
MUT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def get_mut_df(vcf_file, exclude_pat, include_pat, filter=True):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    exclude = []
    skip = np.zeros(3, dtype=int) # missing, no alt, filtered
    data = [[], []]
    idx =[[], []]

    nan_filter = 1
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
            line_reads = np.zeros((sample_no, 4)) # Order: ACGT
            diff_muts = []
            elem_ids = line_cols[8].split(':')

            for s_i, s_rec in enumerate(line_cols[9:]):
                s_rec_elms = s_rec.split(':')
                gt = s_rec_elms[elem_ids.index('GT')].replace('/', '|')
                try:
                    true_gt = s_rec_elms[elem_ids.index('TG')]
                # Missing in Monovar output format
                except ValueError:
                    true_gt = ''
                    line_muts[2, s_i] = np.nan

                # Get true gt
                if true_gt.count('0') == 1:
                    line_muts[2][s_i] = 1
                elif true_gt == '1|1' or true_gt == '2|2' or true_gt == '1|2' \
                        or true_gt == '2|1':
                    line_muts[2][s_i] = 2

                # Set filtered genotypes to nan
                if gt == '.|.':
                    line_muts[0:2, s_i] = np.nan
                    line_reads[s_i] = np.full(4, np.nan)
                    continue
                elif gt == '0|0': # Wildtype
                    pass
                # Check if gt is mutation
                elif '0' in gt: # Heterozygous
                    line_muts[0][s_i] = int(gt[-1])
                    diff_muts.append(int(gt[-1]))
                else: # Homozygous
                    line_muts[1][s_i] = int(gt[-1])
                    diff_muts.append(int(gt[-1]))

                if 'RC' in elem_ids:
                    # Add reads
                    try:
                        reads_i = np.array(
                            s_rec_elms[elem_ids.index('RC')].split(','), dtype=int)
                    except ValueError:
                        reads_i = np.full(4, np.nan)
                    else:
                        if sum(reads_i) == 0:
                            reads_i = np.full(4, np.nan)

                elif 'AD' in elem_ids:
                    try:
                        ref_cnt, alt_cnt = s_rec_elms[elem_ids.index('AD')].split(',')
                    except ValueError:
                        reads_i = np.full(4, np.nan)
                    else:
                        reads_i = np.zeros(4)
                        reads_i[MUT[line_cols[3]]] = int(ref_cnt)
                        reads_i[MUT[line_cols[4]]] = int(alt_cnt)
                        if sum(reads_i) == 0:
                            reads_i = np.full(4, np.nan)

                line_reads[s_i] = reads_i

            # skip lines without any call, neither true or false
            if not any(np.nansum(line_muts, axis=1) > 0):
                skip[1] += 1
                continue
            # Skip lines with >50% missing values (mainly biological data)
            if np.isnan(line_muts[1]).sum() > sample_no * nan_filter:
                skip[0] += 1
                continue

            # Check if called mutations are all the same
            alts, alt_counts = np.unique(diff_muts, return_counts=True)

            # No mutation called
            if alts.size == 0:
                data[0].append(line_muts[0] + line_muts[1])
                alt_idx = MUT[line_cols[4].split(',')[0]] # Take first alt, doesnt matter
            # Two different single mutations: skip
            elif alts.size == 2 and all(alt_counts == 1):
                data[0].append(np.where(np.isnan(line_muts[0]), np.nan, 0))
                line_cols[6] = 'ISA_violation'
            # one single mutation
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
                skip[2] += 1

    muts = pd.DataFrame(data[0], index=idx[1], columns=sample_names)
    true_muts = pd.DataFrame(data[1], index=idx[1], columns=sample_names)

    # Sanity check: remove sample without name
    exclude.append('')
    include_id = [i for i in sample_names if i not in exclude]
    include_idx = [i for i,j in enumerate(sample_names) if j not in exclude]

    if filter:
        muts = muts.loc[idx[0], include_id]
        reads = np.array(reads)[:,include_idx]
        stats = get_stats(muts, true_muts.loc[idx[0], include_id])
    else:
        muts = muts[include_id]
        reads = np.array(reads)
        stats = get_stats(muts, true_muts)

    read_data = (np.array(idx[0]), np.array(ref_gt), reads, np.array(include_id))

    print(f'VCF rows skipped:\n\t{skip[0]: >5} (>f{nan_filter*100}% missing)\n\t'
        f'{skip[1]: >5} (only wt)\n\t{skip[2]: >5} (filtered)')

    return muts, true_muts[include_id], read_data, stats


def get_stats(df, true_df, binary=True, per_cell=False, verbose=False):
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


def show_tree(tree, dendro=False, br_length='mut_no_soft', w_idx=0):
    import copy
    import matplotlib.pyplot as plt
    from matplotlib import cm

    tree = copy.deepcopy(tree)
    tree.ladderize(reverse=False)

    if hasattr(tree.root, 'mut_no_true'):
        true_mut = True
    else:
        true_mut = False
    if hasattr(tree.root, 'weights_norm'):
        weights = True
    else:
        weights = False
    if hasattr(tree.root, 'mut_no_opt'):
        opt_mut = True
    else:
        opt_mut = False

    if true_mut and weights and opt_mut:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}>{c.mut_no_opt:.0f}\n' \
            f'({c.mut_no_true},w{c.weights_norm[w_idx]:.2f})'
    elif true_mut and weights:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}\n' \
            f'({c.mut_no_true},w{c.weights_norm[w_idx]:.2f})'
    elif true_mut and opt_mut:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}>{c.mut_no_opt:.0f}\n' \
            f'({c.mut_no_true})'
    elif weights and opt_mut:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}\n' \
            f'({c.mut_no_true},w{c.weights_norm[w_idx]:.2f})'
    elif true_mut:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}\n({c.mut_no_true})'
    elif weights:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}\n' \
            '(w{c.weights_norm[w_idx]:.2f})'
    elif opt_mut:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}>{c.mut_no_opt:.0f}'
    else:
        br_labels = lambda c: f'{c.mut_no_soft:.1f}'

    try:
        cmap =cm.get_cmap('RdYlBu_r')
        for i in tree.find_clades():
            i.color = [int(j * 255) for j in cmap(i.weights_norm_z[w_idx])[:3]]
    except AttributeError:
        pass

    if br_length == 'lambda':
        tree.root.lambd = 'root'
        br_labels = lambda c: c.lambd
    elif br_length == 'length':
        br_labels = lambda c: c.branch_length
    elif br_length == 'mut_no_soft':
        for i in tree.find_clades():
            i.branch_length = i.mut_no_soft
    elif br_length == 'mut_no_opt':
        for i in tree.find_clades():
            i.branch_length = i.mut_no_opt
    elif br_length == 'mut_no_true':
        for i in tree.find_clades():
            i.branch_length = i.mut_no_true
    else:
        raise RuntimeError(f'Unknown branch length parameter: {br_length}')

    if dendro:
        max_depth = max(tree.depths().values())
        for i in tree.find_clades(terminal=False, order='postorder'):
            int_len = i.name.count('+')
            for child in i.clades:
                child.branch_length = int_len - child.name.count('+')

    plt.rc('font', **{'size': 14})
    plt.rc('lines', **{'linewidth': 3})

    fig = plt.figure(figsize=(10, 20), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, branch_labels=br_labels,
        label_func=lambda c: c.name if c.name and c.name.count('+') == 0 else '',
        subplots_adjust=({'left': 0.01, 'bottom': 0.01, 'right': 0.99, 'top': 0.99}),
        axes=ax
    )


def get_tree(tree_file, paup_exe, samples=[], FN_fix=None, FP_fix=None):
    # Extract data from tree file paths
    try:
        FP = max(float(re.search('WGA0[\.\d,]*-0[\.\d]*-(0[\.\d]*)', tree_file) \
            .group(1)) + 0.01, LAMBDA_MIN) # + 0.01 = Seq. error
        FN = max(float(re.search('WGA(0[\.\d]*)[,\.\d]*?-', tree_file).group(1)),
            LAMBDA_MIN)
    except AttributeError:
        FP = LAMBDA_MIN
        FN = LAMBDA_MIN

    if tree_file.endswith('.raxml.bestTree') \
            or tree_file.endswith('.Mapped.raxml.mutationMapTree') \
            or 'cellphy' in tree_file:
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True,
            br_length=True)

        log_file = tree_file.replace('.mutationMapTree', '.log') \
            .replace('.bestTree', '.log')
        with open(log_file, 'r') as f:
            log = f.read().strip()

        try:
            FP = max(float(re.search('SEQ_ERROR: (0.\d+(e-\d+)?)', log).group(1)),
                LAMBDA_MIN)
            FN = max(float(re.search('ADO_RATE: (0.\d+(e-\d+)?)', log).group(1)),
                LAMBDA_MIN)
        except AttributeError:
            pass

    elif tree_file.endswith('_ml0.newick') or 'scite' in tree_file:
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
            FN = max(LAMBDA_MIN, float(
                re.search('best value for beta:\\\\t(\d.\d+(e-\d+)?)', log_raw) \
                    .group(1)))
            FP = max(LAMBDA_MIN, float(
                re.search('best value for alpha:\\\\t(\d.\d+(e-\d+)?)', log_raw) \
                    .group(1)))
        # For Scite, multiply by two as FN != ADO event if assuming binary data
        FN *= 2
    else:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            br_length=True)

    if FP_fix and FN_fix:
        FP = FP_fix
        FN = FN_fix

    tree = Phylo.read(StringIO(tree_str), 'newick')
    # Prune outgroup (of simulations)
    try:
        outg = [i for i in tree.find_clades() if i.name == 'healthycell'][0]
        tree.prune(outg)
        outg_name = outg.name
    except (IndexError, AttributeError):
        outg_name = None

    # Remove terminals that are not in vcf (discrepancy vcf and newick)
    for node in tree.get_terminals():
        if node.name not in samples:
            tree.prune(node)

    # Initialize node attributes on tree
    for i, node in enumerate(tree.find_clades()):
        node.muts_br = set([])
        node.mut_no_soft = 0
        node.mut_no_hard = 0
        node.name = '+'.join(sorted([j.name for j in node.get_terminals()]))

    return tree, outg_name, FP, FN


def get_tree_reads(tree_file, call_data, paup_exe, w_max, FN_fix=None, FP_fix=None):
    reads = call_data[2]

    tree, outg, FP, FN = get_tree(
        tree_file=tree_file,
        paup_exe=paup_exe,
        samples=reads[3],
        FN_fix=FN_fix,
        FP_fix=FP_fix,
    )
    if outg in reads[3]:
        # Remove outg cell
        cell_idx = np.argwhere(reads[3] != outg).flatten()
        # Remove loci where a mut was only called in the outg
        mut_idx = np.argwhere(
            (call_data[0].drop(outg, axis=1).sum(axis=1) > 0).values).flatten()
        reads = (reads[0][mut_idx], reads[1][mut_idx],
            reads[2][mut_idx][:,cell_idx,:], reads[3][cell_idx])
        # Muts that were only called in outg

    MS = 0
    MS = np.isnan(reads[2].sum(axis=2)).sum(axis=1).sum() \
        / (reads[0].size * reads[3].size)

    M = map_mutations_reads(tree, reads, FP, FN + MS, call_data[1])
    add_br_weights(tree, FP, FN + MS, w_max)

    return tree, FP, FN, M


def get_tree_gt(tree_file, call_data, paup_exe, w_max, FN_fix=None, FP_fix=None):
    muts = call_data[0]

    tree, outg, FP, FN = get_tree(
        tree_file=tree_file,
        paup_exe=paup_exe,
        samples=muts.columns.values,
        FN_fix=FN_fix,
        FP_fix=FP_fix
    )
    FN /= 2

    if outg in muts.columns:
        muts.drop(outg, axis=1, inplace=True)
    # Remove mutations that were only present in outgroup
    muts = muts[muts.sum(axis=1) > 0]

    # Make sure that at FN + MS is max. 0.8
    MS = min(muts.isna().sum().sum() / muts.size, 1 - FN - 0.2)
    M = map_mutations_gt(tree, muts, FP, FN + MS)
    add_br_weights(tree, FP, FN + MS, w_max)

    return tree, FP, FN, M


def _normalize_log_probs(probs, return_normal=True):
    max_i = np.nanargmax(probs)
    try:
        exp_probs = np.exp(probs[np.arange(probs.size) != max_i] \
            - probs[max_i])
    except FloatingPointError:
        exp_probs = np.exp(
            np.clip(probs[np.arange(probs.size) != max_i] - probs[max_i],
                log_LAMBDA_MIN, 0)
        )
    probs_norm = probs - probs[max_i] - np.log1p(np.nansum(exp_probs))

    if return_normal:
        return np.exp(np.clip(probs_norm, log_LAMBDA_MIN, 0))
    else:
        return probs_norm


def _get_S(tree, cells, add_zeros=True):
    n = cells.size
    node_map = {}
    S = pd.DataFrame(data=np.zeros((2 * n - 1, n)), columns=cells)
    # Get leaf nodes of each internal node
    for i, node in enumerate(tree.find_clades()):
        node_map[i] = node
        leaf_nodes = node.name.split('+')
        S.loc[i, leaf_nodes] = 1

    S = S.values
    if add_zeros:
        S = np.vstack([S, np.zeros(n)])
    S_inv = 1 - S
    return S, S_inv, node_map


def map_mutations_gt(tree, muts_in, FP, FN):
    S, S_inv, node_map = _get_S(tree, muts_in.columns)

    idx_map = muts_in.index.values
    nans = np.isnan(muts_in).values
    # Convert ternary hetero/homozygous data to binary mutation data
    muts = np.where(nans, np.nan, muts_in.values.astype(bool).astype(float))
    muts_inv = 1 - muts

    errors = np.log([1 - FN, FP, 1 - FP, FN])
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

        best_nodes = np.argwhere(probs == np.max(probs)).flatten()

        for best_node in best_nodes:
            node_map[best_node].muts_br.add(idx_map[i])
            node_map[best_node].mut_no_hard += 1 / best_nodes.size

    soft_assigned = M.sum(axis=0)
    for i, node in node_map.items():
        node.mut_no_soft = soft_assigned[i]

    root_muts = tree.root.mut_no_soft
    for node in tree.get_terminals():
        node.mut_no_soft_total = root_muts \
            + np.sum([i.mut_no_soft for i in tree.get_path(node)])

    cols = [j.name for i,j in sorted(node_map.items())] + ['FP']
    return pd.DataFrame(M, index=muts_in.index, columns=cols)


def map_mutations_reads(tree, read_data, FP, FN, true_muts):
    from scipy.stats import binom
    # Remove outgroup mutations
    idx, gt, reads, cells = read_data

    S, S_inv, node_map = _get_S(tree, cells)

    gamma = np.log(FN / 2)
    gamma_not = np.log(1 - FN)
    eps = np.log(FP / 3)
    eps_not = np.log(1 - FP)

    eps_het = np.log(2 * FP / 3)
    eps_het_not = np.log(0.5 - FP / 3)

    wt_col = S.shape[0]
    M = np.zeros((reads.shape[0], wt_col))
    skip = np.full(reads.shape[0], False)
    clip_min = -50
    # Get depth and reciproce of depth
    depth = np.sum(reads, axis=2)
    rec_depth = 1 / depth
    # rec_depth = np.ones(depth.size)

    for i in tqdm(range(idx.size)):
        # counts: ref, !ref
        ref = reads[i][:,gt[i][0]]
        ref_not = depth[i] - ref
        # counts: alt, !alt
        alt = reads[i][:,gt[i][1]]
        alt_not = depth[i] - alt
        # counts: ref|alt, !(ref|alt)
        het = ref + alt
        het_not = depth[i] - het
        # Get log prob p(b|0,0) for wt and p(b|1,1) homoyzgous genotypes
        # Clip to avoid that single cells outweights whole prob
        p_count_dist_het = np.clip(
            binom.logpmf(np.nan_to_num(alt), np.nan_to_num(het), 0.5), clip_min, 0)

        p_ref = np.clip((ref * eps_not + ref_not * eps) * rec_depth[i], clip_min, 0)
        p_alt = np.clip((alt * eps_not + alt_not * eps) * rec_depth[i], clip_min, 0)
        # Get log pro for p(b|0,1) for het genotype
        noNAN = ~np.isnan(p_ref)

        p_noADO = np.clip(
            (het * eps_het_not + het_not * eps_het) * rec_depth[i], clip_min, 0)
        p_ADO = np.clip(
            np.logaddexp(gamma + p_ref, gamma + p_alt, where=noNAN), clip_min, 0)

        p_het = np.logaddexp(gamma_not + p_count_dist_het + p_noADO, p_ADO,
            where=noNAN)
        p_mut = np.logaddexp(p_het, p_alt, where=noNAN)

        probs = np.nansum(S_inv * p_ref + S * p_mut, axis=1)
        probs_norm = _normalize_log_probs(probs)

        # Soft assignment
        M[i] = probs_norm
        # Hard assignment
        best_nodes = np.argwhere(probs == np.max(probs)).flatten()
        for best_node in sorted(best_nodes, reverse=True):
            try:
                node_map[best_node].muts_br.add(idx[i])
                node_map[best_node].mut_no_hard += 1 / best_nodes.size
            # Remove muts where MLE is full wt (outside of tree/all zero)
            except KeyError:
                if best_nodes.size == 1:
                    skip[i] = True
                # break

    # Remove muts where MLE contains full wt (outside of tree/all zero)
    print(f'Mutations where ML includes wt: {skip.sum()}')
    skip = np.full(skip.size, False)
    # ------------------------------------------------
    M_filtered = M[~skip]

    soft_assigned = M_filtered.sum(axis=0)
    for i, node in node_map.items():
        node.mut_no_soft = soft_assigned[i]

    cols = [j.name for i,j in sorted(node_map.items())] + ['FP']

    return pd.DataFrame(M_filtered, index=idx[~skip], columns=cols)


def add_br_weights(tree, FP, FN, w_max):
    m = tree.count_terminals()
    n = 2 * m - 1

    # Get successor matrix
    nodes = [i for i in tree.find_clades()]
    leaf_names = sorted([i.name for i in tree.get_terminals()])
    leaf_map = {j: i for i, j in enumerate(leaf_names)}

    S = np.zeros((n, m), dtype=int)
    for i, node in enumerate(nodes):
        cells = [leaf_map[i] for i in node.name.split('+')]
        S[i, cells] = 1

    # Add zero line for wildtype probabilities
    S = np.vstack([S, np.zeros(m)])
    S_inv = 1 - S

    l_TN = np.log(1 - FP)
    l_FN = np.log(FN)
    weights = np.zeros((n, 2), dtype=float)

    for i, y in enumerate(S):
        if i >= n:
            continue

        t = y.sum()
        p_ADO = np.exp(t * l_FN + (m - t) * l_TN)
        p_noADO = 1 - p_ADO

        # weight 0: inv variance
        weights[i, 0] = min(w_max, 1 / (p_noADO * p_ADO))
        # weight 1: odds ratio
        weights[i, 1] = min(w_max, p_noADO / p_ADO)

    # Drop root weight
    root_id = np.argmax(S.sum(axis=1))
    weights = np.delete(weights, root_id, axis=0)
    nodes.remove(tree.root)

    # Normalize to sum to branch number
    weights_norm = weights.shape[0] * weights / weights.sum(axis=0)
    # Normalize for coloring: <= 1 to [0, 0.5], >1 to (0.5, 1]
    weights_norm_z = np.where(weights_norm <= 1, weights_norm / 2,
        0.5 + (weights_norm / (weights_norm.max(axis=0) * 0.5)))

    for i, node in enumerate(nodes):
        node.weights = weights[i]
        node.weights_norm = weights_norm[i]
        node.weights_norm_z = weights_norm_z[i]
    tree.root.weights = np.full(weights.shape[1], -1)
    tree.root.weights_norm = np.full(weights.shape[1], -1)
    tree.root.weights_norm_z = np.ones(weights_norm_z.shape[1])


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
        child0.muts_br_true = child0.true_muts_node \
            .difference(int_node.true_muts_node)
        child0.muts_br_false = child0.muts_br.difference(child0.muts_br_true)
        child0.muts_br_missing = child0.muts_br_true.difference(child0.muts_br)
        child0.mut_no_true = len(child0.muts_br_true)

        child1.muts_br_true = child1.true_muts_node \
            .difference(int_node.true_muts_node)
        child1.muts_br_false = child1.muts_br.difference(child1.muts_br_true)
        child1.muts_br_missing = child1.muts_br_true.difference(child1.muts_br)
        child1.mut_no_true = len(child1.muts_br_true)

    tree.root.muts_br_true = tree.root.true_muts_node
    tree.root.muts_br_false = tree.root.muts_br.difference(tree.root.muts_br_true)
    tree.root.muts_br_missing = tree.root.muts_br_true.difference(tree.root.muts_br)
    tree.root.mut_no_true = len(tree.root.muts_br_true)


def get_sign(x):
    if x < 0:
        return '-'
    else:
        return '+'


def get_model_data(tree, pseudo_mut=0, true_data=False):
    # Lambdas are assigned branchwise from the tips to the root, following the
    #   classical strict molecular clock

    internal_nodes = [i for i in tree.find_clades(terminal=False)]

    leaf_no = tree.count_terminals()
    br_no = 2 * leaf_no - 2

    # Dimensions: 0 soft assignment; 1 hard assignment
    Y = np.zeros((2, br_no), dtype=float)
    init = np.zeros(br_no, dtype=float)
    weights_norm = np.zeros((br_no, tree.root.weights.size), dtype=float)
    constr = np.zeros((leaf_no - 1, br_no), dtype=int)
    constr_cols = []

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
            w0 = node.clades[0].weights_norm[0]
            w1 = node.clades[1].weights_norm[0]
            # Chose traversal with higher weight
            if w0 < w1:
                traverse[br_cells[node.clades[0]]] = 1
                node = node.clades[0]
            else:
                traverse[br_cells[node.clades[1]]] = 1
                node = node.clades[1]
        return traverse

    # Iterate over internal nodes
    node_idx = 0
    br_cells = {}

    def add_node(node, l_idx, neighbor=None):
        br_cells[node] = node_idx
        if true_data:
            Y[1, node_idx] = node.mut_no_true + pseudo_mut
        else:
            Y[1, node_idx] = node.mut_no_hard + pseudo_mut
        Y[0, node_idx] = node.mut_no_soft + pseudo_mut

        weights_norm[node_idx] = node.weights_norm

        if neighbor == None:
            constr[l_idx, node_idx] = 1
            node.lambd = f'+{l_idx}'
        else:
            constr[l_idx] = get_traversal(neighbor, tree) \
                - get_traversal(node, tree)
            node.lambd = ' '.join(sorted([f'{get_sign(j)}{i}' \
                for i, j in enumerate(constr[l_idx]) if j != 0]))

    for l_idx, (_, int_node) in enumerate(sorted_br):
        is_terminal = [i.is_terminal() for i in int_node.clades]

        # Both terminal nodes
        if all(is_terminal):
            # Assign new lambda to terminals
            add_node(int_node.clades[0], l_idx)
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)
            constr_cols.append(int_node.clades[0].name)
            node_idx += 1

            add_node(int_node.clades[1], l_idx)
            constr[l_idx, node_idx] = -1
            constr_cols.append(int_node.clades[1].name)

        # One internal, one terminal node
        elif any(is_terminal):
            # Assign new lambda to terminal
            terminal = int_node.clades[is_terminal.index(True)]
            add_node(terminal, l_idx)
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)
            constr_cols.append(terminal.name)
            node_idx += 1

            # Assign lambda sum to internal
            internal = int_node.clades[is_terminal.index(False)]
            add_node(internal, l_idx, neighbor=terminal)
            constr_cols.append(internal.name)

        # Both internal nodes
        else:
            shorter, longer = sorted(int_node.clades, key=lambda x: x.branch_length)
            # Assign new lambda to shorter branch

            add_node(shorter, l_idx)
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)
            constr_cols.append(shorter.name)
            node_idx += 1

            # Assign lambda sum to longer
            add_node(longer, l_idx, neighbor=shorter)
            constr_cols.append(longer.name)

        init_val = np.dot(constr[l_idx], init)

        # Make sure init is within bounds
        if init_val >= LAMBDA_MIN:
            init[node_idx] = init_val
        else:
            # If not, increase latested added lambda value
            init[node_idx - 1] += -init_val + (0.1)
            init[node_idx] = np.dot(constr[l_idx], init)
        node_idx += 1


    assert np.allclose(constr @ init, 0, atol=LAMBDA_MIN), \
        'Constraints not fulfilled for x_0'
    assert (init >= LAMBDA_MIN).all(), \
        f'Init value smaller than min. distance: {init.min()}'

    # from plotting import plot_weights
    # y_labels = ['ADO', 'Topology', 'Soft']
    # plot_weights(weights_norm, 1, bin_no=weights_norm.shape[0],y_labels=y_labels)

    if true_data:
        Y_out =  Y[1]
    else:
        Y_out = Y[0]

    return Y_out, constr, init, weights_norm, np.array(constr_cols)


def get_LRT_poisson(Y, constr, init, weights=np.array([]), alg='trust-constr'):
    ll_H1 = np.sum((Y * np.log(np.where(Y > 0, Y, 1)) - Y) * weights)
    for i in [10000, 1000, 100, 10, 1, 0.1]:
        scale_fac = (ll_H1 // i) * i
        if scale_fac != 0:
            break

    def fun_opt(l, Y):
        return -np.nansum(
            (Y * np.log(np.where(l > LAMBDA_MIN, l, np.nan)) - l) * weights) \
            / scale_fac

    def fun_jac(l, Y):
        return -(Y / np.clip(l, LAMBDA_MIN, None) - 1) * weights / scale_fac

    def fun_hess(l, Y):
        return np.identity(Y.size) * (Y / np.clip(l, LAMBDA_MIN, None)**2) * weights / scale_fac

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

    ll_H0 = np.nansum((Y * np.log(opt.x) - opt.x) * weights)
    dof = weights.sum() - constr.shape[0]
    LR = -2 * (ll_H0 - ll_H1)

    on_bound = np.sum(opt.x <= LAMBDA_MIN ** 0.5)
    if on_bound > 0:
        dof_diff = np.arange(on_bound + 1)
        p_vals = np.clip(chi2.sf(LR, dof - dof_diff), 1e-100, 1)
        p_val_weights = np.where(
            np.isnan(p_vals), 0, binomCoeff(on_bound, dof_diff))
        p_val = np.average(np.nan_to_num(p_vals), weights=p_val_weights)
    else:
        p_val = chi2.sf(LR, dof)

    devi = ((Y * np.log(np.where(Y > 0, Y, 1)) - Y) * weights) \
        - ((Y * np.log(opt.x) - opt.x) * weights)

    return LR, dof, on_bound, p_val, zip(opt.x, devi)


def run_poisson_tree_test_simulations(vcf_file, tree_file, out_files, paup_exe,
        w_maxs=[999], exclude='', include=''):
    run = os.path.basename(vcf_file).split('.')[1]
    cols = ['-2logLR', 'dof', 'p-value', 'hypothesis', 'weights']

    call_data = get_mut_df(vcf_file, exclude, include)
    # tree, _, _, M = get_tree_reads(tree_file, call_data, paup_exe)

    if not isinstance(w_maxs, list):
        w_maxs = [w_maxs]
    if not isinstance(out_files, list):
        out_files = [out_files]

    for w_max, out_file in zip(w_maxs, out_files):
        tree, FP, FN, M = get_tree_gt(tree_file, call_data, paup_exe, w_max)
        add_true_muts(tree, call_data[1])

        Y, constr, init, weights_norm, constr_cols = get_model_data(tree,
            true_data=False)

        # save_tree(tree, tree_file + 'mapped.newick')

        # weights: 0 = variance, 1 = odds ratio,
        w_idx = 0
        LR, dof, on_bound, p_val, Y_opt = \
            get_LRT_poisson(Y, constr, init, weights_norm[:, w_idx])
        hyp = int(p_val < 0.05)

        header_str = 'run\tSNVs\tTP\tFP\tTN\tFN\tMS\tMS_T\t'
        header_str += '\t'.join([f'{col}_poissonTree_wMax{w_max}' for col in cols])
        model_str = ''

        stats = call_data[3]
        weight_str = ",".join([str(i) for i in weights_norm[:,w_idx].round(3)])

        model_str += f'{run}\t{call_data[0].shape[0]}\t{stats["TP"]}\t{stats["FP"]}\t' \
            f'{stats["TN"]}\t{stats["FN"]}\t{stats["MS"]}\t{stats["MS_T"]}\t' \
            f'{LR:0>5.2f}\t{dof+on_bound}\t{p_val:.2E}\tH{hyp}\t{weight_str}\n'

        with open(out_file, 'w') as f_out:
            f_out.write(f'{header_str}\n{model_str}')

    # add_opt_values_to_tree(tree, Y, Y_opt)
    # show_tree(tree, br_length='length')
    # import pdb; pdb.set_trace()


def add_opt_values_to_tree(tree, Y, Y_opt):
    Y_map = dict(zip(Y, Y_opt))
    for node in tree.find_clades():
        if tree.root == node:
            node.mut_no_opt = -1
        else:
            node.mut_no_opt = Y_map[node.mut_no_soft][0]
            node.deviation_opt = Y_map[node.mut_no_soft][1]


def save_tree(tree, tree_out, w_idx=0):
    int_nodes = 0
    for node in tree.find_clades():
        if node.is_terminal():
            node.name = f'{node.name}[w:{node.weights_norm[w_idx]:.2f}]'
        else:
            node.name = f'node{int_nodes}[w:{node.weights_norm[w_idx]:.2f}]'
            int_nodes += 1
        node.branch_length = node.mut_no_soft

    Phylo.write(tree, tree_out, 'newick',
        format_branch_length='%1.1f', branch_length_only=True)


def run_poisson_tree_test_biological(vcf_file, tree_file, out_file, paup_exe,
        w_max=999, exclude='', include=''):
    call_data = get_mut_df(vcf_file, exclude, include)
    # tree, FP, FN, M = get_tree_reads(tree_file, call_data, paup_exe)
    tree, FP, FN, M = get_tree_gt(tree_file, call_data, paup_exe, w_max)

    Y, constr, init, weights_norm, constr_cols = get_model_data(tree)

    # weights: 0 = variance, 1 = odds ratio
    w_idx = 0
    LR, dof, on_bound, p_val, _ = \
        get_LRT_poisson(Y, constr, init, weights_norm[:,w_idx])
    hyp = int(p_val < 0.05)

    path_strs = vcf_file.split(os.path.sep)
    try:
        clock_dir_no = path_strs.index('ClockTest')
    except ValueError:
        dataset = 'unknown'
        subset = 'unknown'
        filters = 'unknown'
    else:
        subset = path_strs[clock_dir_no + 1]
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

    model_str = f'{dataset}\t{subset}\t{filters}\t{FN:.4f}\t{FP:.4f}\t' \
        f'{LR:0>5.3f}\t{dof}\t{p_val}\tH{hyp}'

    save_tree(tree, tree_file + 'mapped.newick', w_idx)

    header_str = 'dataset\tsubset\tfilters\tFN\tFP\t-2logLR\tdof\tp-value\t' \
        'hypothesis'
    with open(out_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')

    print(model_str)
    # NOTE: has to be here, otherwise subprocess calling script fails
    print('success')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in vcf format')
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, default='', nargs='+',
        help='Output file(s). 1 file per w_max value')
    parser.add_argument('-e', '--exe', type=str,
        default='../../paup4a168_ubuntu64', help='Path to PAUP exe.')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    parser.add_argument('-w', '--w_max', type=float, default=999, nargs='+',
        help='Maximum weight value.')
    parser.add_argument('-b', '--biological_data', action='store_true',
        help='Test true data (instead of simulation data).')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''

        w_maxs = []
        for out_file in snakemake.output:
            file_name = os.path.basename(out_file)
            w_maxs.append(
                float(re.search('wMax(\d+[\.\d]*)\.LRT', file_name).group(1)))
        run_poisson_tree_test_simulations(
            vcf_file=snakemake.input.vcf,
            tree_file=snakemake.input.tree,
            out_files=snakemake.output,
            paup_exe=snakemake.params.paup_exe,
            w_maxs=w_maxs,
            exclude=snakemake.params.exclude,
            include=snakemake.params.include
        )
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = args.vcf + '.poissonTree_LRT.tsv'
        if args.biological_data:
            run_poisson_tree_test_biological(
                vcf_file=args.vcf,
                tree_file=args.tree,
                out_files=args.output,
                paup_exe=args.exe,
                w_maxs=args.w_max,
                exclude=args.exclude,
                include=args.include,
            )
        else:
            run_poisson_tree_test_simulations(
                vcf_file=args.vcf,
                tree_file=args.tree,
                out_files=args.output,
                paup_exe=args.exe,
                w_maxs=args.w_max,
                exclude=args.exclude,
                include=args.include,
            )
