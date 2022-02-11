#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns


FONTSIZE = 16
DPI = 300
RUG_HEIGHT = 0.03
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper')
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
#         'lines.markersize': 6.0,
# })

vis_names = {
    'poissondisp': 'Poisson Dispersion',
    'paup': 'PAUP*',
    'poissontree': 'Poisson Tree',
    'cellcoal': 'True Tree',
    '-': '-',
    'scite': 'Scite Tree',
    'cellphy': 'CellPhy Tree',
}
vis_names_short = {
    'poissondisp': 'Poisson Dispersion',
    'paup': 'PAUP*',
    'poissontree': 'Poisson Tree',
    'cellcoal': 'True Tree',
    '-': '-',
    'scite': 'Scite',
    'cellphy': 'CellPhy',
}

COLORS = [
     # Blue     # Green    # Red      # Orange   # Purple
    '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', # dark
    '#A6CEE3', '#B2DF8A', '#FB9A99', '#FDBF6F', '#CAB2D6', # light
    '#62A3CB', '#72BF5B', '#EF5A5A', '#FE9F37', '#9A77B8', # medium
    '#FFFF99', '#B15928', #ugly
]
colors = {
    'Poisson Dispersion':  '#994EA3', # purple
    'PAUP*': '#4FAF4A', # green
    'PAUP* + True Tree': '#4FAF4A', # green
    'PAUP* + CellPhy Tree': '#177512', # darker green
    'PAUP* + SCITE Tree': '#AAE3A7', #  lighter green
    'cellcoal': '#E41A1A', # red
    'cellphy': '#377DB8', # blue
    'CellPhy': '#377DB8', # blue
    'scite': '#FF7F00', # orange
    'Scite': '#FF7F00', # ligther orange
    '-': '#994EA3' # purple
}
poisson_colors = { # red, blue, orange
    'True Tree': ['#E41A1A', '#377DB8', '#FF7F00'], # normal
    'CellPhy Tree': ['#8C0000', '#094D85', '#9B4D00'], # darker
    'SCITE Tree': ['#F04949', '#84B5DE', '#FFB164'] #brigher
}


HUE_ORDER = ['-', 'PAUP*', 'cellcoal', 'cellphy', 'scite']