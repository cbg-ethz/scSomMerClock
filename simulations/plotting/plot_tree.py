#!/usr/bin/env python3

import argparse
import re

from ete3 import Tree, TreeStyle, NodeStyle, ImgFace

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from  matplotlib.colors import rgb2hex


LAMBDA_MIN = 1e-6
COLORS = [
     # Blue     # Green    # Red      # Orange   # Purple
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', # light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]

TICK_FONTSIZE = 8
LABEL_FONTSIZE = 8

vis_names = {'cellcoal': 'True', 'scite': 'SCITE', 'cellphy': 'CellPhy'}
colors = {
    'True': (0.945, 0.4, 0.627),
    'CellPhy': (0.443, 0.396, 0.776),
    'SCITE': (0.02, 0.604, 0.173)
}


def plot_tree(in_file, out_file=None):

    with open(in_file, 'r') as f:
        tree_str = f.read().strip()

    t = Tree(tree_str, format=1, quoted_node_names=True)

    weight_map = []
    for node in t.traverse():
        weight = float(re.search('.*\[w:(-?\d+\.\d+)\]' ,node.name).group(1))
        if weight == -1:
            node.dist = 0
        # Add 1 psydocount: required for circular plotting
        else:
            node.dist += 1
        node.name = re.search('(.*)\[w:-?\d+\.\d+\]' ,node.name).group(1)
        node.support = weight
        weight_map.append((node, weight))
    # Add edge style for each node
    weights = np.array([j for i,j in weight_map])
    weights_z = np.where(weights <= 1, weights / 2,
        0.5 + 0.5 * weights / weights.max())

    cmap =cm.get_cmap('RdYlBu_r')
    for i, node in enumerate(t.traverse()):
        if weights[i] == -1:
            color_hex = '#000000'
        else:
            color_hex = rgb2hex(cmap(weights_z[i]))

        style = NodeStyle()
        style["size"] = 0 # set internal node size to 0
        style["vt_line_color"] = color_hex
        style["hz_line_color"] = color_hex
        style["vt_line_width"] = 4
        style["hz_line_width"] = 4
        node.img_style = style


    ts = TreeStyle()
    ts.mode = 'r' # c = circular, r = rectangular
    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.show_branch_length = True
    ts.margin_left = 0
    ts.margin_right = 0
    ts.margin_top = 0
    ts.margin_bottom = 0
    ts.root_opening_factor = 0
    ts.show_border = True
    # ts.optimal_scale_level = 'full'
    # safe cbar as image
    cmap_file = in_file + '.cmap.png'
    plot_cmap(cmap, weights[weights >= 0].min(), 1.00, weights.max(), cmap_file)

    ts.legend.add_face(ImgFace(cmap_file, width = 600), column=0)
    ts.legend_position = 4
    # t.show(tree_style=ts)

    if not out_file:
        out_file = in_file + '.PhyloTree.pdf'
    t.render(out_file, tree_style=ts, dpi=300, w=183, units="mm")
    

def plot_cmap(cmap, x_min, x_center, x_max, out_file):
    colors = cmap(np.arange(cmap.N))

    fig, ax = plt.subplots(figsize=(6, 2),
        subplot_kw=dict(xticks=[0, 5, 10], yticks=[]))
    cax = ax.imshow([colors], extent=[0, 10, 0, 1])

    ax.set_xticklabels([f'{x_min:.2f}', f'{x_center:.2f}', f'{x_max:.2f}'])

    fig.subplots_adjust(left=0.05, bottom=0, right=0.95, top=1)
    fig.savefig(out_file, dpi=300)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input newick file.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()    
    plot_tree(args.input, args.output)
