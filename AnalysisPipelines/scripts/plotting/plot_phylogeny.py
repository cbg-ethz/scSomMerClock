#!/usr/bin/env python3

import os
import re

from ete3 import Tree, TreeStyle, NodeStyle, ImgFace, TextFace, CircleFace, RectFace
from ete3.parser.newick import NewickError
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, Normalize, rgb2hex
import numpy as np


def show_tree_full(tree, out_file='', w_idx=0, out_type='pdf', br_labels=False,
            leaf_labels=False, expand_root=True):

    sup_vals = np.unique([i.support for i in tree.iter_descendants()])
    # sup_vals = np.ones(1)

    mut_cutoff = np.inf #max(50,
        #np.percentile([i.dist for i in tree.iter_descendants()], 90))

    try:
        cmap = LinearSegmentedColormap.from_list(
            'blRd', [(0, 0, 0), (1, 0, 0)], N=256)
        dist_fracs = [i.dist_fraction for i in tree.iter_descendants()]
        norm = Normalize(vmin=np.nanmin(dist_fracs), vmax=np.nanmax(dist_fracs))
    except AttributeError:
        cmap =cm.get_cmap('inferno') # RdYlBu_r

    lw = 2
    fsize = 5
    for i, node in enumerate(tree.iter_descendants()):
        if hasattr(node, 'plotted'):
            continue

        style = NodeStyle()

        node.dist += 2
        mut_no = node.dist

        if mut_no > mut_cutoff:
            br_break = TextFace(f'||', fsize=fsize*3)
            br_break.margin_right = 25
            br_break.margin_top = -4
            node.add_face(br_break, column=0, position="float-behind")

            br_break_txt = TextFace(f'{node.dist: <3.0f}', fsize=fsize*2)
            br_break_txt.margin_right = 15
            br_break_txt.margin_bottom = 7
            node.add_face(br_break_txt, column=0, position="branch-top")

            node.dist = mut_cutoff
        else:
            if br_labels:
                node.dist = mut_no
                mut = TextFace(f'{mut_no:.1f}', fsize=fsize)
                mut.margin_right = 2
                node.add_face(mut, column=0, position="branch-top")

        if hasattr(node, 'mut_no_true') and node.mut_no_true >= 0:
            mut = TextFace(f'{node.mut_no_true:.0f} true', fsize=fsize - 1)
            node.add_face(mut, column=0, position="branch-top")

        if False:
            if hasattr(node, 'dist_fraction'):
                color_hex = rgb2hex(cmap(norm(node.dist_fraction)))
                style["vt_line_color"] = color_hex
                style["hz_line_color"] = color_hex
            elif hasattr(node, 'weights_norm_z'):
                if np.abs(node.weights_norm_z[w_idx]) == 1:
                    color_hex = '#000000'
                else:
                    color_hex = rgb2hex(cmap(node.weights_norm_z[w_idx]))
                style["vt_line_color"] = color_hex
                style["hz_line_color"] = color_hex

        # Add support
        if sup_vals.size > 2 and not node.is_leaf() and node.support > 50:
            c1 = CircleFace(8, 'black', )
            c1.margin_top = -4
            c1.margin_right = 0
            c1.margin_left = 0
            c1.margin_bottom = 0
            node.add_face(c1, column=0, position='branch-right')

            c2 = CircleFace(7, 'white')
            c2.margin_top = -15
            c2.margin_right = 1
            c2.margin_left = 1
            c1.margin_bottom = 0
            node.add_face(c2, column=0, position='branch-right')

            supp = TextFace(f'{node.support: >3.0f}', fsize=fsize)
            supp.margin_left = 1
            supp.margin_top = -9
            supp.tight_text = True
            node.add_face(supp, column=0, position="branch-right")

        style["size"] = 0 # set internal node size to 0
        style['vt_line_width'] = lw
        style['hz_line_width'] = lw
        node.img_style = style

        if node.is_leaf():
            if re.match('CRC0907-C\d+', node.name) \
                    or re.match('TPLgr\d+', node.name):
                leaf = RectFace(6, 6, 'black', 'black')
            elif re.match('TDT[NS]+\d+', node.name):
                # leaf = RectFace(3, 3, 'black', 'red')
                ta_pic = os.path.join(
                    os.path.dirname(os.path.realpath(__file__)), 'triangle.png')
                leaf = ImgFace(ta_pic, width=8, height=8)
            else:
                leaf = CircleFace(3, 'black') # RectFace(3, 3, 'black', 'black')
            node.add_face(leaf, column=0, position="branch-right")

            # name = TextFace(f'{node.missing:.2f}', fsize=fsize)
            # name.margin_left = 2
            # name.hz_align = 1
            # node.add_face(name, column=1, position="branch-right")
            if leaf_labels:
                name = TextFace(node.name, fsize=fsize)
                name.margin_left = 2
                name.hz_align = 1
                node.add_face(name, column=0, position="branch-right")


        if br_labels and hasattr(node, 'weights_norm'):
            weight = round(node.weights_norm[w_idx], 1)
            if weight >= 0:
                w_face = TextFace(f'{weight}', fsize=1)
                node.add_face(w_face, column=0, position="branch-bottom")

        if hasattr(node, 'drivers'):
            for i, driver in enumerate(node.drivers):
                if i == 0:
                    space_face = RectFace(0, 0, 'red', 'red')
                    rel_driver_no = len([i for i in node.drivers if i[1]])
                    # space_face.margin_top = 20 * rel_driver_no
                    space_face.opacity = 1
                    node.add_face(space_face, column=0, position="branch-right")
                if driver[2] < 0:
                    continue

                driver_face = TextFace(f'{driver[0]}', fsize=fsize*2)
                driver_face.margin_top = 0
                driver_face.margin_right = 0
                driver_face.hz_align = 1

                if driver[1]:
                    driver_face.background.color = 'LightGreen'
                else:
                    continue
                    driver_face.background.color = 'PeachPuff'
                node.add_face(driver_face, column=0, position="float")
        node.plotted = True


    root = tree.get_tree_root()
    if not expand_root:
        root.children[0].dist = 0
    root.img_style = NodeStyle(size=0)
    root.img_style = style

    ts = TreeStyle()
    ts.mode = 'r' # c = circular, r = rectangular
    # ts.rotation = 90
    ts.allow_face_overlap = False
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False

    ts.min_leaf_separation = 10
    ts.branch_vertical_margin = 3
    ts.margin_left = 20
    ts.margin_right = 5
    ts.margin_top = 5
    ts.margin_bottom = 5
    # ts.optimal_scale_level = 'full' # "mid" | "full"
    if tree.get_farthest_leaf()[1] > 25000:
        ts.scale = None
    else:
        ts.scale = 1

    if out_file:
        if not out_file.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            out_file += f'.{out_type}'

        try:
            tree.render(out_file, tree_style=ts, w=3400, dpi=300)
            print(f'Tree written to: {out_file}')
        except:
            tree.render(out_file, dpi=300, h=2000, units='mm')
            print(f'Simple Tree written to: {out_file}')
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
    
