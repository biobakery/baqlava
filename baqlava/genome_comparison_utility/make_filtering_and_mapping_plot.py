#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from matplotlib.patches import Patch
import sys
from pathlib import Path
from matplotlib import gridspec
import argparse


# parse args:
parser = argparse.ArgumentParser(usage="python make_filtering_and_mapping_plot.py [VGB_###] [marker|protein] [MGX|MTX] [<location of tempfiles>] [<location of script's data files (optional)>]")

parser.add_argument("vgb")
parser.add_argument("datatype", choices=["marker", "protein"])
parser.add_argument("dataset", choices=["MGX", "MTX"])
parser.add_argument("tempfileloc", help="Location of BAQLaVa marker or protein tempfiles")
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

# load in reference files:
data_dir = Path(args.data_dir)
baqref = pd.read_csv(data_dir / "BAQLaVa_metadata_genomes.txt.gz", sep="\t")[['Genome Name','Segment_Group']]
baq_markers = pd.read_csv(data_dir / 'idmap_nucleotide_nofilter.txt', sep="\t", names=['name1','name2','length'])
baq_proteins = pd.read_csv(data_dir / 'idmap_protein_nofilter.txt', sep="\t", names=['name1','name2','length'])

def split_mapname(df, lent):
    df1 = df.copy()
    df1 = df1[df1['length']>= lent]
    name = []
    for i in df1['name1']:
        name.append(i.split("_")[0])
    df1['Genome Name'] = name
    return df1
    
baq_markers = split_mapname(baq_markers, 500)
baq_proteins = split_mapname(baq_proteins, 200)

def gather_mapping_tempfiles(loc, dat):
    
    root = Path(loc)

    if not root.exists():
        raise FileNotFoundError(f"Directory does not exist: {root}")

    if dat == 'Markers':
        files = list(root.rglob("*tempfile_markers.txt"))
    else:
        files = list(root.rglob("*tempfile_proteins.txt"))
        
    # turn into a large dataframe
    retdf = pd.DataFrame({})
    for i in files:
        samplename = i.name.split("_tempfile")[0]
        tempdf = pd.read_csv(str(i), sep="\t")
        tempdf['sample'] = samplename
        retdf = pd.concat([retdf, tempdf])
    return retdf

def compute_average_mapping_markers(df, idmap, ref, VGB):
    df1 = pd.merge(idmap.copy(), ref.copy(), on='Genome Name', how='inner')
    df1 = df1[df1['Segment_Group']==VGB].drop_duplicates()
    # for later: get total markerized length per genome:
    ML = df1.copy()[['Genome Name','length']].groupby('Genome Name', as_index=False).sum().rename(columns={'length':'genome_markerized_length'})

    # get all markers expected to see x all samples (so we can include zeros) 
    samplelist_df = pd.DataFrame({"sample": list(set(m1['sample']))})
    expanded_df1 = df1.merge(samplelist_df, how="cross")[['Genome Name','name1','sample','length']]
    
    # now add in the mapping info we found:
    df2 = pd.merge(expanded_df1, df.copy(), left_on=['name1','sample'], right_on=['Marker','sample'], how='left')[['Genome Name','name1','sample','RPK','length']].fillna(0)
    
    # this version average over **markers as boolean** (what percent of samples mapped this marker?)
    df3 = df2.copy()
    df3['RPK'] = df3['RPK'].astype('bool').astype('int')
    df3 = df3[['Genome Name','name1','RPK']].groupby(['Genome Name','name1'], as_index=False).mean()

    # this version average over **genomes as boolean** (what percent of samples mapped this genome?)
    # to do this we want to calculate total mapped length, which baqlava uses as a threshold
    df4 = df2.copy().query("RPK>0")[['Genome Name','sample','length']].groupby(['Genome Name','sample'], as_index=False).sum().rename(columns={'length':'mapped_marker_length'})
    df5 = pd.merge(df4, ML, on='Genome Name', how='left')
    df5['percent_markerized_length_mapped'] = df5['mapped_marker_length'] / df5['genome_markerized_length']
    df5["genome_mapped"] = (df5["percent_markerized_length_mapped"] >= 0.5).astype(int)
    df5 = df5[['Genome Name','genome_mapped']].groupby('Genome Name', as_index=False).mean().rename(columns={'genome_mapped':'percent_samples_mapped'})
    return df3, df5

def compute_average_mapping_proteins(df, idmap, ref, VGB):
    df1 = pd.merge(idmap.copy(), ref.copy(), on='Genome Name', how='inner')
    df1 = df1[df1['Segment_Group']==VGB].drop_duplicates()
    # for later: get total proteome length per genome:
    ML = df1.copy()[['Genome Name','length']].groupby('Genome Name', as_index=False).sum().rename(columns={'length':'genome_proteome_length'})

    # get all proteins expected to see x all samples (so we can include zeros) 
    samplelist_df = pd.DataFrame({"sample": list(set(m1['sample']))})
    expanded_df1 = df1.merge(samplelist_df, how="cross")[['Genome Name','name1','sample','length']]
    
    # now add in the mapping info we found:
    df2 = pd.merge(expanded_df1, df.copy(), left_on=['name1','sample'], right_on=['Protein','sample'], how='left')[['Genome Name','name1','sample','RPK','length']].fillna(0)
    
    # this version average over **proteins as boolean** (what percent of samples mapped this protein?)
    df3 = df2.copy()
    df3['RPK'] = df3['RPK'].astype('bool').astype('int')
    df3 = df3[['Genome Name','name1','RPK']].groupby(['Genome Name','name1'], as_index=False).mean()

    # this version average over **genomes as boolean** (what percent of samples mapped this genome?)
    # to do this we want to calculate percent of ORFs mapped, which baqlava uses as a threshold
    # and also impose the minimum proteome length threshold
    df4 = df2.copy()
    df4['RPK'] = df4['RPK'].astype('bool').astype('int')
    df4 = df4.copy()[['Genome Name','sample','RPK']].groupby(['Genome Name','sample'], as_index=False).mean().rename(columns={'RPK':'percent_proteome_mapped'})
    df4['mapped_percent_threshold'] = df4['percent_proteome_mapped']>=0.5
    df4['mapped_percent_threshold'] = df4['mapped_percent_threshold'].astype('int')
    # add in proteome length:
    df5 = pd.merge(df4.copy(), ML.copy(), on='Genome Name', how='left')
    df5['mapped_length_threshold'] = df5['genome_proteome_length']>=2500
    df5['mapped_length_threshold'] = df5['mapped_length_threshold'].astype('int')
    # get the genomes that pass both  thresholds:
    df6 = df5.copy()
    df6['pass_thresholds'] = df6['mapped_length_threshold'] * df6['mapped_percent_threshold']
    df6 = df6[['Genome Name','pass_thresholds']].groupby('Genome Name', as_index=False).mean().rename(columns={'pass_thresholds':'percent_samples_mapped'})
    return df3, df6

# get mapping information:
m1 = gather_mapping_tempfiles(args.tempfileloc, datatype)

if datatype == 'Markers':
    avg_markers, avg_genomes = compute_average_mapping_markers(m1, baq_markers, baqref, VGB)
elif datatype == 'Proteins':
    avg_markers, avg_genomes = compute_average_mapping_proteins(m1, baq_proteins, baqref, VGB)

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

# final merge of mapping df with plot information:
plotdat3 = pd.merge(plotdat1.copy(), avg_genomes.copy().set_index('Genome Name'), on='Genome Name', how='left').fillna(0)[['percent_samples_mapped']]


# Make the plot:
fig_height = len(plotdat1) * 0.15

fig = plt.figure(figsize=(8, fig_height))

gs = gridspec.GridSpec(
    2, 2,
    width_ratios=[4, 1],                 # left panel wider
    height_ratios=[3, len(plotdat1)],
    wspace=0.15,
    hspace=0.02
)

# Left column
ax_top = fig.add_subplot(gs[0, 0])
ax_bottom = fig.add_subplot(gs[1, 0])

# Bottom-right panel
ax_right = fig.add_subplot(gs[1, 1])


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
    bbox_to_anchor=(1.03, 1.5),
    handlelength=1.5,
    handleheight=1.5,
    labelspacing=1.2)

# -------------------------
# BOTTOM RIGHT MARKER MAPPING
# -------------------------

# Number of genomes / rows
nrows = len(plotdat1)

# y positions matching heatmap row centers
ypos = np.arange(nrows) + 0.5

ax_right.barh(
    ypos,
    plotdat3["percent_samples_mapped"].values,
    height=1.0,
    edgecolor='white',
    linewidth=1,
    color='#8aa6bf'
)

# Make ax_right fill the same vertical span as ax_bottom heatmap
ax_right.set_ylim(ax_bottom.get_ylim())

# Optional cleanup
ax_right.set_xlim(0, 1)
ax_right.set_yticks([])
ax_right.set_xlabel("Percent of samples\nmeeting mapping threshold")

ax_right.spines['top'].set_visible(False)
ax_right.spines['right'].set_visible(False)


fig.savefig(f"{VGB}_{datatype}_{args.dataset}_mapping.png", bbox_inches="tight", dpi=300)
