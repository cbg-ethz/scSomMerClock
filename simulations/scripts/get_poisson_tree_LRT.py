#!/usr/bin/env python3

import os
import re
from io import StringIO

import numpy as np
import pandas as pd
from scipy.stats.distributions import chi2

from utils import change_newick_tree_root

from newick import read, loads
from Bio import Phylo


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

    excl_leafs = []
    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_terminals():
        cols_flag = leaf_node.name == muts.columns

        if any(cols_flag):
            leaf_s = muts.loc[:, cols_flag].squeeze()
            others_s = muts.loc[:, ~cols_flag].sum(axis=1)
            no_muts = ((leaf_s == 1) & (others_s == 0)).sum()
            leaf_node.branch_length = no_muts
        else:
            excl_leafs.append(leaf_node.name)

    # Remove excluded leaf nodes
    for excl_leaf in excl_leafs:
        tree.collapse(excl_leaf)
        
    # Add number of mutations to internal nodes
    for int_node in tree.get_nonterminals():
        leaf_desc = [i.name for i in int_node.get_terminals()]
        cols_flag = muts.columns.isin(leaf_desc)
        leaf_no = len(leaf_desc)
        
        int_s = muts.loc[:, cols_flag].sum(axis=1)
        others_s = muts.loc[:, ~cols_flag].sum(axis=1)

        no_muts = ((int_s == leaf_no) & (others_s == 0)).sum()
        int_node.branch_length = no_muts
        int_node.name = '+'.join(muts.columns[cols_flag])

    tree.ladderize(reverse=False)
    Phylo.write(tree, "example.newick", "newick")

    total_muts = muts.shape[0]
    cell_no = tree.count_terminals()
    br_no = 2 * cell_no - 2
    br_indi_no = cell_no - 1

    tree_mat = np.zeros((br_no, br_indi_no))
    y = []
    br_cells = {}

    n_idx = 0 # node index
    l_idx = 1 # lambda index
    
    sorted_br = sorted(
        [(tree.distance(i), i) for i in tree.find_clades(terminal=False)],
        key=lambda x: x[0], reverse=True
    )

    for l_idx, (mut_no, int_node) in enumerate(sorted_br, 1):
        for child in int_node.clades:
            y.append(child.branch_length)
            br_cells[child] = n_idx
            tree_mat[n_idx, :l_idx] = 1
            if not child.is_terminal():
                tree_mat[n_idx] -= tree_mat[br_cells[child.clades[0]]]
                if not all(tree_mat[br_cells[child.clades[0]]] == tree_mat[br_cells[child.clades[1]]]):
                    import pdb; pdb.set_trace()
            n_idx += 1
        
    import pdb; pdb.set_trace()



    tree_mat = np.zeros((br_no, br_indi_no))
    level = 0
    idx = 0
    curr_level_nodes = tree.descendants
    next_lvl_nodes = []
    while level < br_indi_no:
        for curr_node in curr_level_nodes:
            if curr_node.is_leaf:
                br_cells.append(curr_node.name)
                tree_mat[idx, level:] = 1
                idx += 1
            else:
                next_lvl_nodes = curr_node.descendants
                node_samples = '+'.join(
                    sorted([i.name for i in curr_node.get_leaves()]))
                br_cells.append(node_samples)
                tree_mat[idx, level] = 1
                idx += 1
        curr_level_nodes = next_lvl_nodes
        level += 1
        import pdb; pdb.set_trace()


def test_poisson(vcf_file, tree_file, out_file, paup_exe, exclude='', include='',
            alpha=0.05):

    run = os.path.basename(vcf_file).split('.')[1]
    muts = get_mut_matrix(vcf_file, exclude, include)
    tree = get_tree_dict(tree_file, muts, paup_exe)
    import pdb; pdb.set_trace()
    mean_muts = muts.mean()


    LR = 2 * np.sum(muts * np.log(muts / mean_muts))
    # LR2 =  np.sum((muts - mean_muts)**2) / mean_muts
    dof = len(muts) - 1
    p_val = chi2.sf(LR, dof)
    if p_val < alpha:
        hyp = 'H1'
    else:
        hyp = 'H0'

    out_str = f'{run}\t{LR:0>5.2f}\t{dof}\t{p_val:.2E}\t{hyp}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('run\t-2logLR\tdof\tp-value\thypothesis\n')
        f_out.write(out_str)


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