#!/usr/bin/env python3

from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

DEF_ADO = [0, 0.05, 0.1, 0.2, 0.4, 0.6]
FONTSIZE = 16
DPI = 300
RUG_HEIGHT = 0.03
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'xtick.major.size': 2, 'ytick.major.size': 2, 'lines.linewidth': 2,
        'axes.axisbelow': True})
#     rc={'font.size': FONTSIZE,
#         'axes.labelsize': 'medium',
#         'axes.titlesize': 'large',
#         'xtick.labelsize': 'medium',
#         'ytick.labelsize': 'medium',
#         'legend.fontsize': 'medium',
#         'legend.title_fontsize': 'large',
#         'axes.labelticksize': 50,
#         'lines.linewidth': 1,
#         'xtick.major.size':  6,
#         'ytick.major.size':  6,
#         'lines.markersize': 6.0}
# })


METHODS = {
    'PAUP*': 'PAUP* (True)',
    'poissonDisp': 'Poisson dispersion',
    'cellcoal': 'Poisson tree (True)',
    'cellphy': 'Poisson tree (CellPhy)',
    'scite': 'Poisson tree (Scite)',
    'mobster': 'mobster',
    'neutrality': r'1/$f$ test'
}

COLORS = {
    'PAUP*': '#f4a259',
    'poissonDisp': '#6ec15d',
    'cellcoal': '#e06c78',
    'cellphy': '#5874dc',
    'scite': '#bb4430',
    'neutrality': '#FFE652',
    'mobster': '#2E851C', # mobster
    1: '#2E851C', # mobster
    0: '#384e78', # mobster
}

HUE_ORDER = ['cellcoal', 'cellphy', 'poissonDisp', 'PAUP*', 'scite']
HIST_DEFAULT = {
    'alpha': 0.75,
    'fill': True,
    'binwidth': 0.05,
    'binrange': (0, 1),
    'element': 'bars',
    'stat': 'probability',
    'kde': False,
    'common_norm': False,
    'fill': True,
    'multiple': 'layer',
    'palette': COLORS,
    'hue_order': HUE_ORDER,
    'legend': False,
}
LINE_STYLE ={
    10: (0, (1, 4, 1, 4)), # . . . .
    30: (0, (2, 3, 1, 3)), # - . - .
    50: (0, (2, 3, 2, 3)), # -- . -- .
    70: (0, (4, 2, 2, 2)),
    90: (0, (4, 1, 4, 1)), # -- -- -- --
    100: (0, (1, 0, 0, 0)), # solid
}
MARGINS = {
    'left': 0.15,
    'right': 0.9,
    'top': 0.8,
    'bottom': 0.3,
    'wspace': 0.25,
    'hspace': 0.75,
}

bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=1)


def get_subplots(row_no, col_no, scale=1):
    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no * 1.5 * scale, row_no * 2 * scale))
    axes = np.reshape(axes, (row_no, col_no))
    return fig, axes


def plot_fig(fig=None, out_file=''):
    # fig.tight_layout()
    plt.subplots_adjust(**MARGINS)
    if out_file:
        if not out_file.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            out_file += '.png'
        fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()



def add_rugs(data, offset, ax, color, height=None):
    if not height:
        height = RUG_HEIGHT

    segs = np.stack((np.c_[data, data],
            np.c_[np.zeros_like(data) + 1 + height *2 * offset,
                    np.zeros_like(data) + 1 + height *2 * (offset + 1)]),
        axis=-1)
    lc = LineCollection(segs, transform=ax.get_xaxis_transform(),
        clip_on=False, color=color, linewidth=0.06)
    ax.add_collection(lc)


def add_col_header(ax, title):
    ax.annotate(title, xy=(0.5, 1.1), xytext=(0, 5), xycoords='axes fraction',
        textcoords='offset points', ha='center', va='baseline',
        annotation_clip=False)


def add_row_header(ax, title):
    ax2 = ax.twinx()
    ax2.set_ylabel(title)
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    for tick in  ax.yaxis.majorTicks:
        tick.tick1line.set_markersize(0)