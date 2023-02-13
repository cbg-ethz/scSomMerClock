#!/usr/bin/env python3

import os
import re
from io import StringIO

import numpy as np
# np.seterr(all='raise')
import pandas as pd

from utils import change_newick_tree_root

from ete3 import Tree

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


def get_tree_dict(tree_file, muts, paup_exe):
    if 'cellphy_dir' in tree_file:
        _, tree_str = change_newick_tree_root(tree_file, paup_exe, root=True)
    elif 'trees_dir' in tree_file:
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False)
    elif 'scite_dir' in tree_file:
        samples = [f'tumcell{i:0>4d}' for i in range(1, muts.shape[1], 1)] \
            + ['healthycell']
        tree_str, _ = change_newick_tree_root(tree_file, paup_exe, root=False,
            sample_names=samples)

    # With BioPython package
    try:
        tree = Tree(tree_str, format=2)
    except NewickError:
        tree = Tree(tree_str, format=1)

    # Add number of mutations to terminal nodes
    for leaf_node in tree.get_leaves():
        leaf_muts = muts[leaf_node.name]
        others_muts = muts.drop(leaf_node.name, axis=1).sum(axis=1)
        no_muts = ((leaf_muts == 1) & (others_muts == 0)).sum()
        leaf_node.branch_length = no_muts
        
    # Add number of mutations to internal nodes
    for int_node in tree.iter_descendants():
        if int_node.is_leaf():
            continue
        leaf_desc = [i.name for i in int_node.get_leaves()]
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


def show_pvals(p_vals):
    import matplotlib.pyplot as plt
    
    _plot_pvals(p_vals)
    plt.show()
    plt.close()


def get_mut_distribution(vcf_files, tree_files, out_file, paup_exe, exclude='', include=''):
    Ys = []
    for i, vcf_file in enumerate(vcf_files):

        run = os.path.basename(vcf_file).split('.')[1]
        muts = get_mut_matrix(vcf_file, exclude, include)
        tree = get_tree_dict(tree_files[i], muts, paup_exe)
        Y, X_H0, labels = get_model_data(tree, muts)
        Ys.extend(Y)

    vals, cnts = np.unique(Ys, return_counts=True)
    with open(out_file, 'w') as f:
        f.write('Value\tCounts\n')
        for i, val in enumerate(vals):
            f.write(f'{val}\t{cnts[i]}\n')

    # from plotting import _plot_muts
    # _plot_muts(Ys)
    # plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str, nargs='+', help='Input files')
    parser.add_argument('tree', type=str, nargs='+', help='Tree files in newick format')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
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
        get_mut_distribution(snakemake.input.vcf, snakemake.input.tree, 
            snakemake.output[0], snakemake.params.paup_exe,
            snakemake.params.exclude, snakemake.params.include)
    else:
        import argparse
        args = parse_args()
        get_mut_distribution(args.vcf, args.tree, args.output, args.exe,
            args.exclude, args.include)