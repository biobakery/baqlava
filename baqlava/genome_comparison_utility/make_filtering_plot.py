#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from matplotlib.patches import Patch
import sys
import argparse
from pathlib import Path

# parse args:
parser = argparse.ArgumentParser(usage="python make_plasmid_filtering_plot.py [VGB_###] [marker|protein] [MGX|MTX] [<location of datafiles (optional)>]")

parser.add_argument("vgb")
parser.add_argument("datatype", choices=["marker", "protein"])
parser.add_argument("dataset", choices=["MGX", "MTX"])
parser.add_argument(
    "data_dir",
    nargs="?",
    default=Path(__file__).resolve().parent,
    help="Location of data files (defaults to script directory)"
)

args = parser.parse_args()

# set whether we use marker or protein annotation:
datatype = "Markers" if args.datatype == "marker" else "Proteins"

VGB = args.vgb
VGB_f = VGB.split("_")[-1]

# get just information for VGB of interest:
filepath = args.data_dir / f"{args.datatype}_plasmid_filtering_{args.dataset}.txt.gz"
data = pd.read_csv(filepath, sep="\t")
plotdat = data.copy()
plotdat = plotdat[plotdat['Segment_Group']==VGB]

# Get main marker/protein data:
plotdat1 = plotdat.copy()[['Genome Name','Filter Level: None','Filter Level: Default','Filter Level: Conservative','Filter Level: V1_full_VGBs']].drop_duplicates()
plotdat1['order'] = plotdat1[['Filter Level: None','Filter Level: Default','Filter Level: Conservative','Filter Level: V1_full_VGBs']].sum(axis=1)
plotdat1 = plotdat1.sort_values(by=['order','Genome Name'], ascending=[False, True]).set_index('Genome Name')[['Filter Level: None','Filter Level: V1_full_VGBs','Filter Level: Default','Filter Level: Conservative']]

# Get just VGB-level bar:
plotdat2 = plotdat.copy()[['Segment_Group','Filter Level: None','Filter Level: Default','Filter Level: Conservative','Filter Level: V1_full_VGBs']].groupby('Segment_Group', as_index=False).max()
plotdat2 = plotdat2[['Segment_Group','Filter Level: None','Filter Level: V1_full_VGBs','Filter Level: Default','Filter Level: Conservative']].set_index('Segment_Group')

# Make the plot:
fig_height = len(plotdat1) * 0.15

# Create stacked layout
fig, (ax_top, ax_bottom) = plt.subplots(2,1,
                                        figsize=(6, fig_height),
                                        gridspec_kw={'height_ratios': [3, len(plotdat1)], 'hspace': 0.02})

# Colors
soft_red = '#e57373'
soft_green = '#66bb6a'
cmap = ListedColormap([soft_red, soft_green])

# -------------------------
# TOP SINGLE-ROW HEATMAP
# -------------------------

sns.heatmap(data = plotdat2,
            cmap = cmap,
            linewidths = 0.005,
            cbar = False,
            ax = ax_top,
            vmin = 0, vmax = 1)

# Remove top x labels
ax_top.xaxis.tick_top()
ax_top.set_xticklabels(
    ['None', 'BAQLaVa V1\nVGB-Level', 'Default', 'Conservative'],
    rotation=0,
    ha='center')
ax_top.xaxis.set_label_position('top')

ax_top.set_yticklabels(ax_top.get_yticklabels(),rotation=0, size=12)
ax_top.set_ylabel("")

ax_top.set_title("Filtering Level Applied", size=15, pad=15)

# -------------------------
# BOTTOM MAIN HEATMAP
# -------------------------

sns.heatmap(data = plotdat1, 
            cmap = cmap,
            linewidths = 0.005,
            cbar = False,
            ax = ax_bottom,
            vmin = 0, vmax = 1)

# Force all labels to appear
nrows = len(plotdat1)
ytick_size = max(5, min(10, 1000 / nrows))
ax_bottom.set_yticks(np.arange(len(plotdat1.index)) + 0.5)
ax_bottom.set_yticklabels(plotdat1.index, fontsize=ytick_size)

ax_bottom.set_ylabel(f"VGB {VGB_f} Genomes", size=15)
ax_bottom.set_xticks([])

legend_elements = [
    Patch(facecolor=soft_red, edgecolor='black', label=f"{datatype} are not present\nat this filtering level"),
    Patch(facecolor=soft_green, edgecolor='black', label=f"{datatype} are present\nat this filtering level")]

ax_top.legend(
    handles=legend_elements,
    loc='center left',
    bbox_to_anchor=(1.03, -1),
    handlelength=1.5,
    handleheight=1.5,
    labelspacing=1.2)

fig.savefig(f"{VGB}_{datatype}_{args.dataset}_plasmid_filtering.png", bbox_inches="tight", dpi=300)
