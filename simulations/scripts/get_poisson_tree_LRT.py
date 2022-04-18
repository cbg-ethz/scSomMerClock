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

from ete3 import Tree


LAMBDA_MIN = 1e-6
log_LAMBDA_MIN = -100
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

            pos = (int(line_cols[0]), int(line_cols[1]))

            idx[1].append(pos)
            if 'PASS' in line_cols[6] or line_cols[6] == 's50' \
                    or line_cols[6] == 'singleton':
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
    else:
        muts = muts[include_id]
        reads = np.array(reads)

    read_data = (np.array(idx[0]), np.array(ref_gt), reads, np.array(include_id))

    return muts, true_muts[include_id], read_data


def show_tree(tree, out_file='', w_idx=0, out_type='pdf', br_labels=True):
    from ete3 import TreeStyle, NodeStyle, ImgFace, TextFace
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import rgb2hex

    mut_cutoff = max(50,
        np.percentile([i.dist for i in tree.iter_descendants()], 90))

    cmap =cm.get_cmap('RdYlBu_r')
    lw = 3
    fsize = 4
    for i, node in enumerate(tree.iter_descendants()):
        style = NodeStyle()

        mut_no = node.dist
        if mut_no > mut_cutoff:
            node.dist = mut_cutoff
            style['hz_line_type'] = 1
            mut = TextFace(f'{mut_no:.1f}', fsize=fsize)
            mut.margin_right = 5
            mut.vt_align = 1
            node.add_face(mut, column=0, position="branch-top")
        else:
            if br_labels:
                node.dist = mut_no
                mut = TextFace(f'{mut_no:.1f}', fsize=fsize)
                mut.margin_right = 2
                node.add_face(mut, column=0, position="branch-top")

        if hasattr(node, 'mut_no_true'):
            mut = TextFace(f'{node.mut_no_true:.0f} true', fsize=fsize - 1)
            node.add_face(mut, column=0, position="branch-top")

        if node.weights_norm_z[w_idx] == -1:
            color_hex = '#000000'
        else:
            color_hex = rgb2hex(cmap(node.weights_norm_z[w_idx]))

        style["size"] = 0 # set internal node size to 0
        style["vt_line_color"] = color_hex
        style["hz_line_color"] = color_hex
        style["vt_line_width"] = lw
        style["hz_line_width"] = lw
        node.img_style = style

        if node.is_leaf():
            name = TextFace(node.name, fsize=fsize)
            name.margin_left = 2
            name.hz_align = 1
            node.add_face(name, column=0, position="branch-right")

        if br_labels:
            weight = round(node.weights_norm[w_idx], 1)
            if weight >= 0:
                w_face = TextFace(f'{weight}', fsize=1)
                node.add_face(w_face, column=0, position="branch-bottom")

        for i, driver in enumerate(node.drivers):
            if driver[2] < 0:
                continue

            driver_face = TextFace(f'{driver[0]} ({driver[2]* 100:.0f}%)', fsize=3)
            driver_face.hz_align = 1

            if driver[1]:
                driver_face.background.color = 'LightGreen'
            else:
                continue
                driver_face.background.color = 'PeachPuff'
            node.add_face(driver_face, column=0, position="branch-bottom")

    root = tree.get_tree_root()
    root.children[0].dist = 0
    root.img_style = NodeStyle(size=0)

    ts = TreeStyle()
    ts.mode = 'r' # c = circular, r = rectangular
    # ts.arc_span = 180
    # ts.root_opening_factor = 0.1
    # ts.rotation = 90
    ts.allow_face_overlap = True
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.branch_vertical_margin = 5
    ts.min_leaf_separation = 2
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.margin_bottom = 10
    # ts.optimal_scale_level = 'full' # "mid" | "full"
    ts.scale = 1

    leaves = tree.get_leaves()
    leaf_root_dist = np.zeros(len(leaves))
    for i, node in enumerate(leaves):
        while node:
           leaf_root_dist[i] += node.dist
           node = node.up
    h = (max(leaf_root_dist) - min(leaf_root_dist)) / 400
    # ts.scale = 1 / max(leaf_root_dist) * 5

    if out_file:
        # ts.show_border = True
        # Plot cMap to file
        # cmap_file = os.path.splitext(out_file)[0] + '.cmap.png'
        # colors = cmap(np.arange(cmap.N))
        # fig, ax = plt.subplots(figsize=(5, 1),
        #     subplot_kw=dict(xticks=[0, 5, 10], yticks=[]))
        # cax = ax.imshow([colors], extent=[0, 10, 0, 1])
        # fig.subplots_adjust(left=0.1, bottom=0, right=0.9, top=1.5)
        # fig.savefig(cmap_file, dpi=300)

        # ts.legend.add_face(ImgFace(cmap_file, height=40), column=0)
        # ts.legend_position = 4

        if not out_file.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            out_file += f'.{out_type}'

        try:
            tree.render(out_file, tree_style=ts, w=2400,  dpi=300)
            print(f'Tree written to: {out_file}')
        except:
            tree.render(out_file, dpi=300, h=h, units='mm')
            print(f'Simple Tree written to: {out_file}')
    else:
        tree.show(tree_style=ts)


def read_tree(tree_file, samples=[]):
    with open(tree_file, 'r') as f:
        tree_raw = f.read().strip()

    if tree_file.endswith('.raxml.bestTree') \
            or tree_file.endswith('.Mapped.raxml.mutationMapTree') \
            or 'cellphy' in tree_file:

        tree_raw = re.sub('\[\d+\]', '', tree_raw)
        log_file = tree_file.replace('.mutationMapTree', '.log') \
            .replace('.bestTree', '.log')

        if log_file == tree_file:
            log_file = tree_file.replace('newick', 'log')

        with open(log_file, 'r') as f:
            log = f.read().strip()

        FP = float(re.search('SEQ_ERROR: (0.\d+(e-\d+)?)', log).group(1))
        FN = float(re.search('ADO_RATE: (0.\d+(e-\d+)?)', log).group(1))
    elif tree_file.endswith('_ml0.newick') or 'scite' in tree_file:
        sem_count = tree_raw.count(';')
        if sem_count == 0:
            tree_raw += ';'
        elif sem_count > 1:
            tree_raw = tree_raw.split(';')[0].strip() + ';'

        nodes = [int(i) for i in \
            re.findall('(?<=[\(\),])\d+(?=[,\)\(;])', tree_raw)]
        int_node_idx = 1
        for i, s_i in enumerate(sorted(nodes)):
            if i < samples.size:
                pat = f'(?<=[\(\),]){s_i}(?=[,\)\(;)])'
                repl = f'{samples[i]}:1.0'
            else:
                pat = f'(?<=[\(\),]){s_i}(?=[,\)\(;)])'
                repl = f'Node{int_node_idx}:1.0'
                int_node_idx += 1
            tree_raw = re.sub(pat, repl, tree_raw)
        tree_raw = tree_raw[tree_raw.index('('):]

        # Get ML error rates from log file
        log_file =  tree_file.replace('_ml0.newick', '.log')
        if log_file == tree_file:
            log_file =  tree_file.replace('.newick', '.log')

        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                log_raw = f.read()
            FN = float(
                re.search('best value for beta:\\\\t(\d.\d+(e-\d+)?)', log_raw) \
                    .group(1))
            FP = float(
                re.search('best value for alpha:\\\\t(\d.\d+(e-\d+)?)', log_raw) \
                    .group(1))
        else:
            error_file = tree_file.replace('_ml0.newick', '.errors.csv')
            with open(error_file, 'r') as f:
                errors_raw = f.read().strip().split('\n')
            FP, FN = [float(i) for i in errors_raw[1].split(',')]

        # For Scite, multiply by two as FN != ADO event if assuming binary data
        FN *= 2
    else:
        tree_raw = re.sub('cell(?=\d+)', 'tumcell', tree_raw)
        tree_raw = re.sub('outgcell', 'healthycell', tree_raw)

        # Extract data from tree file paths
        try:
            FP = float(re.search('WGA0[\.\d,]*-0[\.\d]*-(0[\.\d]*)', tree_file) \
                .group(1)) + 0.01 # + 0.01 = Seq. error

            FN = float(re.search('WGA(0[\.\d]*)[,\.\d]*?-', tree_file).group(1))
        except AttributeError:
            FP = LAMBDA_MIN
            FN = LAMBDA_MIN

    outg_name = 'healthycell'
    tree = Tree(tree_raw, format=1)
    outg_node = tree&outg_name
    tree.set_outgroup(outg_node)
    outg_node.delete()

    # Remove terminals that are not in vcf (discrepancy vcf and newick; subsampling)
    for node in tree.get_leaves():
        if node.name not in samples:
            node.delete()

    # Initialize node attributes on tree
    weight_no = 2
    for i, node in enumerate(tree.iter_descendants()):
        node.add_features(
            in_length=node.dist,
            muts_br=[[], [], []],
            muts_br_true=set([]),
            muts_br_false=set([]),
            muts_br_missing=set([]),
            muts_node_true=set([]),
            mut_no_soft=0,
            mut_no_true=-1,
            name='+'.join(sorted(node.get_leaf_names())),
            weights=np.zeros(weight_no),
            weights_norm=np.zeros(weight_no),
            weights_norm_z=np.zeros(weight_no),
            drivers=[]
        )

    return tree, outg_name, FP, FN


def get_gt_tree(tree_file, call_data, w_max, FN_fix=None, FP_fix=None):
    muts = call_data[0]

    tree, outg, FP, FN = read_tree(tree_file, muts.columns.values)
    if FP_fix and FN_fix:
        FP =  max(FP_fix, LAMBDA_MIN)
        FN =  max(FN_fix, LAMBDA_MIN)
    else:
        FP = max(FP, LAMBDA_MIN)
        FN /= 2
        FN = max(FN, LAMBDA_MIN)

    # Remove mutations that were only present in outgroup
    if outg in muts.columns:
        muts_red = muts.drop(outg, axis=1)
    muts_red = muts_red[muts_red.sum(axis=1) > 0]

    # if not (tree_file.endswith('.raxml.bestTree') \
    #         or tree_file.endswith('.Mapped.raxml.mutationMapTree') \
    #         or tree_file.endswith('_ml0.newick') \
    #         or 'scite' in tree_file or 'cellphy' in tree_file):
    #     true_red = call_data[1].drop(outg, axis=1)

    #     FN = max(((true_red == 1) & (muts_red == 0)).mean().mean(), LAMBDA_MIN)
    #     FP = max(((true_red == 0) & (muts_red == 1)).mean().mean(), LAMBDA_MIN)

    # Make sure that at FN + MS is max. 0.8
    MS = min(muts_red.isna().mean().mean(), 1 - FN - 0.2)
    MS_i = np.clip(muts_red.isna().mean(), None, 1 - FN - 0.2)
    M = map_mutations_gt(tree, muts_red, FP, FN + MS_i)

    add_br_weights(tree, FP, FN, MS_i, w_max)

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


def map_mutations_gt(tree, muts_in, FP, FN):
    n = muts_in.shape[1]
    node_map = {}
    S = pd.DataFrame(data=np.zeros((2 * n - 1, n)), columns=muts_in.columns)
    # Get leaf nodes of each internal node
    for i, node in enumerate(tree.iter_descendants()):
        node_map[i] = node
        leaf_nodes = node.name.split('+')
        S.loc[i, leaf_nodes] = 1

    S = S.values
    S = np.vstack([S, np.zeros(n)])
    S_inv = 1 - S

    idx_map = muts_in.index.values
    nans = np.isnan(muts_in).values
    # Convert ternary hetero/homozygous data to binary mutation data
    muts = np.where(nans, np.nan, muts_in.values.astype(bool).astype(float))
    muts_inv = 1 - muts

    # -- Assume average missing value for all cells --
    # # errors = np.log([1 - FN.mean(), FP, 1 - FP, FN.mean()])
    # TP_m_all = S * np.log(1 - FN.mean())
    # FP_m_all = S_inv * np.log(FP)
    # TN_m_all = S_inv * np.log(1 - FP)
    # FN_m_all = S * np.log(FN.mean())

    # -- Taking different missing value per cell into account --
    FN_i = FN[muts_in.columns].values

    TP_m = S * np.log(1 - FN_i)
    FP_m = S_inv * np.log(FP)
    TN_m = S_inv * np.log(1 - FP)
    FN_m = S * np.log(FN_i)

    M = np.zeros((muts_in.shape[0], S.shape[0]))
    for i, mut in enumerate(muts):
        # probs = np.stack([
        #     np.nansum(TP_m_all * mut, axis=1), # TP
        #     np.nansum(FP_m_all * mut, axis=1), # FP
        #     np.nansum(TN_m_all * muts_inv[i], axis=1), # TN
        #     np.nansum(FN_m_all * muts_inv[i], axis=1), # FN,
        # ]).sum(axis=0)

        probs = np.stack([
            np.nansum(TP_m * mut, axis=1), # TP
            np.nansum(FP_m * mut, axis=1), # FP
            np.nansum(TN_m * muts_inv[i], axis=1), # TN
            np.nansum(FN_m * muts_inv[i], axis=1), # FN,
        ]).sum(axis=0)

        probs_norm = _normalize_log_probs(probs)
        M[i] = probs_norm

        best_nodes = np.argwhere(probs == np.max(probs)).flatten()
        max_prob = probs_norm.max()
        unique =  best_nodes.size == 1

        for best_node in best_nodes:
            node_map[best_node].muts_br[0].append(idx_map[i])
            node_map[best_node].muts_br[1].append(max_prob)
            node_map[best_node].muts_br[2].append(unique)

    soft_assigned = M.sum(axis=0)
    for i, node in node_map.items():
        node.mut_no_soft = soft_assigned[i]
        node.dist = soft_assigned[i]

    cols = [j.name for i,j in sorted(node_map.items())] + ['FP']
    return pd.DataFrame(M, index=muts_in.index, columns=cols)


def add_br_weights(tree, FP, FN, MS_i, w_max):
    m = len(tree)
    n = 2 * m - 1

    # Get successor matrix
    nodes = [i for i in tree.iter_descendants()]
    leaf_names = sorted([i.name for i in tree.get_leaves()])
    leaf_map = {j: i for i, j in enumerate(leaf_names)}

    S = np.zeros((n, m), dtype=int)
    for i, node in enumerate(nodes):
        cells = [leaf_map[i] for i in node.name.split('+')]
        S[i, cells] = 1

    l_TN = np.log(1 - FP)
    weights = np.zeros((n, 2), dtype=float)

    # # -- Assume average missing value for all cells --
    # l_FN = np.log(FN + MS_i[leaf_names].mean())
    # t = S.sum(axis=1)
    # p_ADO_all = np.exp(t * l_FN.mean() + (m - t) * l_TN)
    # p_noADO_all = 1 - p_ADO_all

    # -- Taking different missing value per cell into account --
    l_FN = np.log(FN)
    MS_vals = MS_i[leaf_names].values

    p_ADO = np.zeros(n, dtype=float)
    for i, br in enumerate(S):
        p_br = log_LAMBDA_MIN

        l_mut = br.sum()
        mut_coeffs = np.log(binomCoeff(l_mut, np.arange(l_mut + 1)))
        l_mut_MS_avg = np.log(MS_vals[br].mean())

        l_no_mut = m - l_mut
        no_mut_coeffs = np.log(binomCoeff(l_no_mut, np.arange(l_no_mut + 1)))
        l_no_mut_MS_avg = np.log(MS_vals[1 - br].mean())

        for FN_no, mut_coeff in enumerate(mut_coeffs):
            MS_mut_no = l_mut - FN_no

            tbl_prob = log_LAMBDA_MIN
            for TN_no, no_mut_coeff in enumerate(no_mut_coeffs):
                MS_no_mut_no = l_no_mut - TN_no

                tbl_prob_new = no_mut_coeff \
                    + FN_no * l_FN + MS_mut_no * l_mut_MS_avg \
                    + TN_no * l_TN + MS_no_mut_no * l_no_mut_MS_avg
                tbl_prob = np.logaddexp(tbl_prob, tbl_prob_new)


            p_br = tbl_prob = np.logaddexp(p_br, mut_coeff + tbl_prob_new)
        p_ADO[i] = np.exp(p_br)
    p_noADO = 1 - p_ADO

    # weight 0: inv variance
    weights[:,0] = np.clip(1 / (p_noADO * p_ADO), None, w_max)
    # weight 1: odds ratio
    weights[:,1] = np.clip(p_noADO / p_ADO, None, w_max)

    # Drop root weight
    root_id = np.argwhere(S.sum(axis=1) == m).flatten()[0]
    weights = np.delete(weights, root_id, axis=0)
    nodes.remove(nodes[root_id])

    # Normalize to sum to branch number
    weights_norm = weights.shape[0] * weights / weights.sum(axis=0)
    # Normalize for coloring: <= 1 to [0, 0.5], >1 to (0.5, 1]
    weights_norm_z = np.where(weights_norm <= 1, weights_norm / 2,
        0.5 + (weights_norm / (weights_norm.max(axis=0) * 0.5)))

    for i, node in enumerate(nodes):
        node.weights = weights[i]
        node.weights_norm = weights_norm[i]
        node.weights_norm_z = weights_norm_z[i]
    root = tree.get_tree_root().children[0]
    root.weights = np.full(weights.shape[1], -1)
    root.weights_norm = np.full(weights.shape[1], -1)
    root.weights_norm_z = np.ones(weights_norm_z.shape[1])


def add_true_muts(tree, df_true):
    for node in tree.iter_descendants('postorder'):
        if node.is_leaf():
            node.muts_node_true = set(df_true.index[df_true[node.name] > 0])
            continue

        child0 = node.children[0]
        child1 = node.children[1]
        node.muts_node_true = child0.muts_node_true \
            .intersection(child1.muts_node_true)
        # child0.false_muts = node.muts_node_true.difference(node.true_muts)
        child0.muts_br_true = child0.muts_node_true \
            .difference(node.muts_node_true)
        child0.muts_br_false = set(child0.muts_br[0]) \
            .difference(child0.muts_br_true)
        child0.muts_br_missing = child0.muts_br_true \
            .difference(set(child0.muts_br[0]))
        child0.mut_no_true = len(child0.muts_br_true)

        child1.muts_br_true = child1.muts_node_true \
            .difference(node.muts_node_true)
        child1.muts_br_false = set(child1.muts_br[0]) \
            .difference(child1.muts_br_true)
        child1.muts_br_missing = child1.muts_br_true \
            .difference(set(child1.muts_br[0]))
        child1.mut_no_true = len(child1.muts_br_true)

    root = tree.get_tree_root().children[0]
    root.muts_br_true = root.muts_node_true
    root.muts_br_false = set(root.muts_br[0]).difference(root.muts_br_true)
    root.muts_br_missing = root.muts_br_true.difference(set(root.muts_br[0]))
    root.mut_no_true = len(root.muts_br_true)


def add_opt_muts(tree, Y, Y_opt):
    Y_map = dict(zip(Y, Y_opt))
    for node in tree.iter_descendants():
        if len(tree) == len(node):
            node.mut_no_opt = -1
        else:
            node.mut_no_opt = Y_map[node.mut_no_soft][0]


def get_model_data(tree, pseudo_mut=0, true_data=False):
    # Lambdas are assigned branchwise from the tips to the root, following the
    #   classical strict molecular clock

    leaf_no = len(tree)
    br_no = 2 * leaf_no - 2

    # Dimensions: 0 soft assignment; 1 hard assignment
    Y = np.zeros((2, br_no), dtype=float)
    init = np.zeros(br_no, dtype=float)
    weights_norm = np.zeros((br_no, tree.get_tree_root().children[0].weights.size))
    constr = np.zeros((leaf_no - 1, br_no), dtype=int)
    constr_cols = []
    br_cells = {}

    def get_traversal(start_node, tree):
        traverse = np.zeros(br_no, dtype=int)

        traverse[br_cells[start_node]] = 1
        node = start_node
        while True:
            # Exit case: terminal node
            if len(node.children) == 0:
                break
            w0 = node.children[0].weights_norm[0]
            w1 = node.children[1].weights_norm[0]
            # Chose traversal with higher weight
            if w0 < w1:
                traverse[br_cells[node.children[0]]] = 1
                node = node.children[0]
            else:
                traverse[br_cells[node.children[1]]] = 1
                node = node.children[1]
        return traverse

    def add_node(node, l_idx, neighbor=None):
        br_cells[node] = node_idx

        Y[0, node_idx] = node.mut_no_soft + pseudo_mut
        Y[1, node_idx] = node.mut_no_true + pseudo_mut

        weights_norm[node_idx] = node.weights_norm

        if neighbor == None:
            constr[l_idx, node_idx] = 1
        else:
            constr[l_idx] = get_traversal(neighbor, tree) \
                - get_traversal(node, tree)

    # Iterate over internal nodes
    node_idx = 0
    int_nodes = [i for i in tree.iter_descendants('postorder') if not i.is_leaf()]
    for l_idx, int_node in enumerate(int_nodes):
        is_terminal = [i.is_leaf() for i in int_node.children]

        # Both terminal nodes
        if all(is_terminal):
            # Assign new lambda to terminals
            add_node(int_node.children[0], l_idx)
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)
            constr_cols.append(int_node.children[0].name)
            node_idx += 1

            add_node(int_node.children[1], l_idx)
            constr[l_idx, node_idx] = -1
            constr_cols.append(int_node.children[1].name)

        # One internal, one terminal node
        elif any(is_terminal):
            # Assign new lambda to terminal
            terminal = int_node.children[is_terminal.index(True)]
            add_node(terminal, l_idx)
            init[node_idx] = max(Y[0, node_idx] + pseudo_mut, LAMBDA_MIN)
            constr_cols.append(terminal.name)
            node_idx += 1

            # Assign lambda sum to internal
            internal = int_node.children[is_terminal.index(False)]
            add_node(internal, l_idx, neighbor=terminal)
            constr_cols.append(internal.name)

        # Both internal nodes
        else:
            shorter, longer = sorted(int_node.children, key=lambda x: x.mut_no_soft)
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
    dof = round(weights.sum() - constr.shape[0])
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


def run_poisson_tree_test_simulations(vcf_file, tree_file, out_file, w_maxs,
            exclude='', include='', FN_fix=None, FP_fix=None):
    run = os.path.basename(vcf_file).split('.')[1]
    cols_all = ['FN', 'FP']
    cols_wMax = ['-2logLR', 'dof', 'p-value', 'hypothesis', 'weights']

    call_data = get_mut_df(vcf_file, exclude, include)

    if not isinstance(w_maxs, list):
        w_maxs = [w_maxs]

    header_str = 'run\t'
    model_str = f'{run}\t'
    for i, w_max in enumerate(sorted(w_maxs)):
        tree, FP, FN, M = get_gt_tree(tree_file, call_data, w_max, FN_fix, FP_fix)
        add_true_muts(tree, call_data[1])

        Y, constr, init, weights_norm, constr_cols = get_model_data(tree)

        # weights: 0 = variance, 1 = odds ratio,
        w_idx = 0
        LR, dof, on_bound, p_val, Y_opt = \
            get_LRT_poisson(Y, constr, init, weights_norm[:, w_idx])

        if i == 0:
            header_str += f'FN\tFP\t'
            model_str += f'{FN:.4f}\t{FP:.6f}\t'

        weight_str = ",".join([str(j) for j in weights_norm[:,w_idx].round(3)])
        header_str += '\t'.join([f'{j}_poissonTree_wMax{w_max}' for j in cols_wMax])
        model_str += f'{LR:0>5.2f}\t{dof+on_bound}\t{p_val:.2E}\t' \
            f'H{int(p_val < 0.05)}\t{weight_str}'

    with open(out_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')
    # add_opt_muts(tree, Y, Y_opt)
    # show_tree(tree)
    # import pdb; pdb.set_trace()


def run_poisson_tree_test_biological(vcf_file, tree_file, out_file,
        w_maxs=[1000], plot_only=False, driver_file='', exclude='', include=''):
    path_strs = vcf_file.split(os.path.sep)
    try:
        clock_dir_no = path_strs.index('ClockTest')
    except ValueError:
        dataset = path_strs[-1].split('_')[0]
        if dataset == 'S21':
            dataset = '_'.join(path_strs[-1].split('_')[:2])
        subset = 'unknown'
        filters = 'unknown'
    else:
        subset = path_strs[clock_dir_no + 1]
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

    if driver_file:
        drivers = pd.read_csv(driver_file, sep='\t')
        drivers.set_index('SYMBOL', inplace=True)

        maf_file = os.path.join(os.path.dirname(vcf_file), f'{dataset}.maf')
        annot = pd.read_csv(maf_file, sep='\t', skiprows=[0], usecols=
            ['SYMBOL', 'Chromosome', 'Start_Position', 'End_Position',
                 'Variant_Classification'])
        annot.set_index(['Chromosome', 'Start_Position'], inplace=True)
        # annot.dropna(how='all',axis=1, inplace=True)
        # If statistics required on the annotations, this is the place
        annot.dropna(inplace=True)
        annot = annot[~annot['Variant_Classification'].str.contains('Intron|Silent')]

    header_str = 'dataset\tsubset\tfilters\tFN\tFP'
    model_str = ''
    w_cols = ['-2logLR', 'dof', 'p-value', 'hypothesis']

    call_data = get_mut_df(vcf_file, exclude, include)
    for w_max in w_maxs:
        tree, FP, FN, M = get_gt_tree(tree_file, call_data, w_max)
        Y, constr, init, weights_norm, constr_cols = get_model_data(tree)

        # weights: 0 = variance, 1 = odds ratio
        w_idx = 0
        LR, dof, on_bound, p_val, _ = \
            get_LRT_poisson(Y, constr, init, weights_norm[:,w_idx])
        hyp = int(p_val < 0.05)

        header_str += '\t' + '\t'.join([f'{i}_{w_max:.0f}' for i in w_cols])
        model_str += f'\t{LR:0>5.3f}\t{dof}\t{p_val}\tH{hyp}'

        if plot_only:
            annotate_drivers(dataset, tree, drivers, annot)
            tree_fig = out_file + f'_w{w_max:.0f}_mapped'
            show_tree(tree, tree_fig, w_idx)
            exit()

    model_str = f'{dataset}\t{subset}\t{filters}\t{FN:.4f}\t{FP:.4f}' + model_str

    summary_file = re.sub('(?<=_)w\d+_', '', out_file)
    with open(summary_file, 'w') as f_out:
        f_out.write(f'{header_str}\n{model_str}')

    # NOTE: has to be here, otherwise subprocess calling script fails
    print('success')


def annotate_drivers(dataset, tree, drivers, annot):
    CANCER_TYPE = {
        'H65': 'MM',
        'Li55': 'BLCA',
        'Lo-P1': None, 'Lo-P2': None, 'Lo-P3': None,
        'Ni8': 'LUAD',
        'S21_P1': 'PRAD', 'S21_P2': 'PRAD', 'S21': 'PRAD',
        'W32': 'BRCA', 'W55': 'BRCA',
        'Wu61': 'COREAD', 'Wu63': 'COREAD',
        'X25': 'RCCC',
    }

    print(f'{dataset}: {CANCER_TYPE[dataset]}')

    for node in tree.iter_descendants():
        for i, (chr, pos) in enumerate(node.muts_br[0]):
            try:
                SNV = annot.loc[chr, pos]
            except TypeError:
                continue

            gene = SNV['SYMBOL']
            try:
                driver = drivers.loc[gene]
            except KeyError:
                continue

            # Skip drivers mapped to more than 1 branch
            if not node.muts_br[2][i] and node.muts_br[2][i] < 0.33:
                continue

            driver_str = f'{SNV["Variant_Classification"]: <17} in '\
                f'driver gene {gene: <15}'
            if isinstance(driver, pd.Series):
                if driver['CANCER_TYPE'] == CANCER_TYPE[dataset]:
                    c_type = True
                    driver_str = '  * ' + driver_str
                    driver_str += f' ({driver["ROLE"]})'
                else:
                    c_type = False
                    driver_str = '    ' + driver_str
                    driver_str += f' (cancer type: {driver["CANCER_TYPE"]})'
            else:
                c_types = driver['CANCER_TYPE'].unique()
                if CANCER_TYPE[dataset] in c_types:
                    c_type = True
                    role = driver[driver['CANCER_TYPE'] == CANCER_TYPE[dataset]] \
                        ['ROLE'].unique()
                    driver_str = '  * ' + driver_str
                    driver_str += f' ({",".join(role)})'
                else:
                    c_type = False
                    driver_str = '    ' + driver_str
                    driver_str += f' (cancer types: {",".join(c_types)})'
            print(driver_str)

            node.drivers.append((gene, c_type, node.muts_br[1][i]))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, help='SNP file in vcf format')
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    parser.add_argument('-excl', '--exclude', type=str, default='',
        help='Regex pattern for samples to exclude from LRT test,')
    parser.add_argument('-incl', '--include', type=str, default='',
        help='Regex pattern for samples to include from LRT test,')
    parser.add_argument('-w', '--w_max', type=float, nargs='+',
        default=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
        help='Maximum weight value. Defaut = 100, 200, ..., 1000')
    parser.add_argument('-b', '--biological_data', action='store_true',
        help='Test true data (instead of simulation data).')
    parser.add_argument('-p', '--plotting', action='store_true',
        help='Plot phylogeny of biological data only.')
    parser.add_argument('-d', '--drivers', type=str, default='',
        help='Path to IntOGen driver genes for annotation of real data sets.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.params.exclude == None:
            snakemake.params.exclude = ''
        if snakemake.params.include == None:
            snakemake.params.include = ''

        run_poisson_tree_test_simulations(
            vcf_file=snakemake.input.vcf,
            tree_file=snakemake.input.tree,
            out_file=snakemake.output[0],
            w_maxs=snakemake.params.w_maxs,
            exclude=snakemake.params.exclude,
            include=snakemake.params.include,
            FN_fix=snakemake.params.FN,
            FP_fix=snakemake.params.FP,
        )
    else:
        import argparse
        args = parse_args()
        if not args.output:
            args.output = [args.vcf + '.poissonTree_LRT.tsv']
        if args.biological_data:
            run_poisson_tree_test_biological(
                vcf_file=args.vcf,
                tree_file=args.tree,
                out_file=args.output,
                w_maxs=args.w_max,
                plot_only=args.plotting,
                driver_file=args.drivers,
                exclude=args.exclude,
                include=args.include,
            )
        else:
            run_poisson_tree_test_simulations(
                vcf_file=args.vcf,
                tree_file=args.tree,
                out_file=args.output,
                w_maxs=args.w_max,
                exclude=args.exclude,
                include=args.include,
            )
