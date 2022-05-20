#!/usr/bin/env python3

import os
import re

from ete3 import Tree, TreeStyle, NodeStyle, ImgFace, TextFace, CircleFace
from ete3.parser.newick import NewickError
import numpy as np


def show_tree(tree, out_file='', out_type='pdf'):
    sup_vals = np.unique([i.support for i in tree.iter_descendants()])

    br_length_all = [i.dist for i in tree.iter_descendants()]
    br_length_cutoff = max(50, np.percentile(br_length_all, 90))

    lw = 1
    fsize = 4
    for i, node in enumerate(tree.iter_descendants()):
        style = NodeStyle()

        if node.dist > br_length_cutoff:
            style['hz_line_type'] = 1
            br_length = TextFace(f'{node.dist:.2f}', fsize=fsize)
            br_length.margin_right = 5
            br_length.vt_align = 1
            node.add_face(br_length, column=0, position="branch-top")

        # Add support
        if sup_vals.size > 1 and not node.is_leaf():
            c1 = CircleFace(4, 'black', )
            c1.margin_top = -3
            c1.margin_right = -2
            c1.margin_left = 0
            c1.margin_bottom = 0
            node.add_face(c1, column=0, position='branch-right')

            c2 = CircleFace(3.5, 'white')
            c2.margin_top = -7.5
            c2.margin_right = -2
            c2.margin_left = 0.5
            c1.margin_bottom = 0
            node.add_face(c2, column=0, position='branch-right')

            supp = TextFace(f'{node.support: >3.0f}', fsize=fsize-2)
            supp.margin_left = 1
            supp.margin_top = -5.5
            supp.tight_text = True
            node.add_face(supp, column=0, position="branch-right")

        style["size"] = 0 # set internal node size to 0
        style['vt_line_width'] = lw
        style['hz_line_width'] = lw
        node.img_style = style

        if node.is_leaf():
            name = TextFace(node.name, fsize=fsize)
            name.margin_left = 2
            name.hz_align = 1
            node.add_face(name, column=0, position="branch-right")

    root = tree.get_tree_root()
    style = NodeStyle()
    style["size"] = 0 # set internal node size to 0
    root.img_style = style

    ts = TreeStyle()
    ts.mode = 'r'
    ts.allow_face_overlap = True
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.branch_vertical_margin = 2
    ts.min_leaf_separation = 2
    ts.margin_left = 10
    ts.margin_right = 10
    ts.margin_top = 10
    ts.margin_bottom = 10
    # ts.scale = 1

    if out_file:
        if not out_file.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            out_file += f'.{out_type}'

        tree.render(out_file, tree_style=ts, dpi=300, w=1600, units='px')
        print(f'Tree written to: {out_file}')
    else:
        tree.show(tree_style=ts)


def read_tree(tree_file, samples=[]):
    with open(tree_file, 'r') as f:
        tree_raw = f.read().strip()

    cellphy_ends = ('.raxml.bestTree', 'raxml.mutationMapTree',
        'raxml.supportFBP', 'raxml.supportTBE')
    if tree_file.endswith(cellphy_ends) or 'cellphy' in tree_file:
        tree_raw = re.sub('\[\d+\]', '', tree_raw)
        if tree_raw.count('tumcell') == 0:
            tree_raw = re.sub('cell(?=\d+)', 'tumcell', tree_raw)
        if tree_raw.count('healthycell') == 0:
            tree_raw = re.sub('outgcell', 'healthycell', tree_raw)

    elif tree_file.endswith('_ml0.newick') or 'scite' in tree_file:
        sem_count = tree_raw.count(';')
        if sem_count == 0:
            tree_raw += ';'
        elif sem_count > 1:
            tree_raw = tree_raw.split(';')[0].strip() + ';'
        tree_raw = tree_raw[tree_raw.index('('):]
    else:
        if tree_raw.count('tumcell') == 0:
            tree_raw = re.sub('cell(?=\d+)', 'tumcell', tree_raw)
        if tree_raw.count('healthycell') == 0:
            tree_raw = re.sub('outgcell', 'healthycell', tree_raw)


    try:
        tree = Tree(tree_raw, format=2)
    except NewickError:
        tree = Tree(tree_raw, format=1)

    if tree_file.endswith('_ml0.newick') or 'scite' in tree_file:
        outg_name = str(max([int(i.name) for i in tree.get_leaves()]))
    else:
        outg_name = 'healthycell'

    outg_node = tree&outg_name
    anc_node = outg_node.get_ancestors()[0]
    if anc_node.is_root():
        outg_node.support = 100
    else:
        outg_node.support = anc_node.support
    tree.set_outgroup(outg_node)
    outg_node.delete()

    return tree


def plot_tree(tree_file, out_file):
    tree = read_tree(tree_file)
    show_tree(tree, out_file)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tree', type=str, help='Tree file in newick format')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    import argparse
    args = parse_args()
    if os.path.isdir(args.tree):
        for file in os.listdir(args.tree):
            if not file.endswith('newick') or 'scite' in file:
                continue
            in_file = os.path.join(args.tree, file)
            out_file = f'{in_file}.phylogeny.png'
            if False: #os.path.exists(out_file):
                print(f'!WARNING : output exists, skipping: {out_file}')
                continue
            print(file)
            plot_tree(in_file, out_file)
    else:
        plot_tree(args.tree, args.output)
    
