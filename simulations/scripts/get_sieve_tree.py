#!/usr/bin/env python3

import re
from utils import change_newick_tree_root


def get_sieve_tree(in_file, out_file, cells, paup_exe=''):
    if 'cellphy_dir' in in_file:
        _, tree = change_newick_tree_root(in_file, paup_exe, root=True,
            br_length=True)
    elif 'trees_dir' in in_file:
        tree, _ = change_newick_tree_root(in_file, paup_exe, root=False,
            br_length=True)
    elif 'scite_dir' in in_file:
        samples = [f'tumcell{i:0>4d}' for i in range(1, cells + 1, 1)]
        tree, _ = change_newick_tree_root(in_file, paup_exe, root=False,
            sample_names=samples, br_length=True)

    tree_cut = re.sub('(,healthycell:0.\d+)|(healthycell:0.\d+,)', '', tree)

    with open(out_file, 'w') as f_out:
        f_out.write(tree_cut)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-c', '--cells', type=int, help='Number of tumor cells.')
    parser.add_argument('-e', '--exe', type=str, help='Path to PAUP exe.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        get_sieve_tree(snakemake.input[0], snakemake.output[0],
            snakemake.params.cells, snakemake.params.paup_exe)
    else:
        import argparse
        args = parse_args()
        get_sieve_tree(args.input, args.output, args.cells, args.exe)