import pandas as pd
import os
import sys
import numpy as np
import scipy
from scipy import stats
from Bio import SeqIO


# Config Parsers
# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser

try:
    import ConfigParser as configparser
except ImportError:
    import configparser

config = configparser.ConfigParser()

install_folder=os.path.dirname(os.path.realpath(__file__))
config_file=os.path.join(install_folder,os.path.pardir,"configs/baqlava.cfg")
config.read(config_file)

nucleotide_reference_file = os.path.abspath(config.get('utility','nucleotide_reference'))
protein_reference_file = os.path.abspath(config.get('utility','protein_reference'))
VGB_taxonomy_file = os.path.abspath(config.get('utility','VGB_taxonomy'))

nucleotide_reference = pd.read_csv(nucleotide_reference_file, sep="\t", low_memory=False)
protein_reference = pd.read_csv(protein_reference_file, sep="\t", low_memory=False)
VGB_taxonomy = pd.read_csv(VGB_taxonomy_file, sep="\t", index_col='Unnamed: 0', low_memory=False)

def format_file_and_base(in_genefams):
    # get column of interest (basename plus RPK)
    base = ".".join(os.path.split( in_genefams )[1].split(".")[:-1]).replace('_genefamilies', '_Abundance-RPKs')
    # remove duplicate abundance lines (since everything is unclassified)
    df1 = pd.read_csv(in_genefams, sep="\t").copy()
    df2 = df1[~df1['# Gene Family'].str.contains("\|")]
    df3 = df2[df2['# Gene Family'] != 'UNMAPPED']
    return base, df3


def get_read_length(file):
    counter = 0
    lens = []
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            counter += 1
            if counter <=10000:
                lens.append(len(record.seq))
            else:
                break
    return sum(lens) / len(lens)

def process_baqlava_nucleotide(base, format_df, ref, readlen, taxref):
    newbase = "_".join(base.split("_nucleotide_"))
    
    df1 = format_df.copy()
    df2 = pd.merge(df1, ref.copy()[['VGB','marker','len','cluster_member', 'segment_group']], left_on='# Gene Family', right_on='marker', how='right')
    df2[base] = df2[base].fillna(0)    
    df2['reads_mapped'] = df2[base] * (df2['len']/1000) 
    df3 = df2[df2['len']>400]
    df3 = df3.copy()
    df3['marker_len_adj'] = df3['len'] - (readlen-1) 

    # get VGBs/segment groups with >75% markerized length mapped to:
    # (for VGBs, 75% of one cluster representative genome's markerized length)
    # (for segment groups, 75% of the 'core' VGB markerized length) 
    # VGBs first:
    df31 = df3.copy()[['VGB','segment_group','marker','cluster_member','len',base]].drop_duplicates()
    df32 = df31.copy().groupby(['VGB','cluster_member'], as_index=False).sum()[["VGB","cluster_member","len"]]
    df33 = df31.copy()[df31[base]!=0].groupby(['VGB','cluster_member'], as_index=False).sum()[["VGB","cluster_member","len"]]
    df34 = pd.merge(df32, df33, on=['cluster_member','VGB'], how='outer').fillna(0)
    df34['p_marker_len_mapped'] = df34['len_y'] / df34['len_x']
    df35 = df34[df34["p_marker_len_mapped"]>=0.75][['VGB']].drop_duplicates()
    df36 = pd.merge(df35, ref.copy()[['VGB','segment_group']].drop_duplicates(), on='VGB', how='left')
    df37 = df36[~df36["segment_group"].str.contains("segment")][['segment_group']]
    # now segments: 
    df31s = pd.merge(df31.copy(), ref.copy()[['VGB','conserved_VGB']].drop_duplicates(), on='VGB', how='left').fillna("not_segmented").query("conserved_VGB!='not_segmented'")
    df32s = df31s.copy().groupby(["segment_group","conserved_VGB"], as_index=False).sum().query("conserved_VGB=='conserved'")[["segment_group","len"]]
    df33s = df31s.copy()[df31s[base]!=0]
    df34s = df33s.groupby(["segment_group","conserved_VGB"], as_index=False).sum().query("conserved_VGB=='conserved'")[["segment_group","len"]]
    df35s = pd.merge(df32s, df34s, on='segment_group', how='left').fillna(0)
    df35s['p_conserved_marker_len_mapped'] = df35s['len_y'] / df35s['len_x']
    df36s = df35s[df35s["p_conserved_marker_len_mapped"]>=0.75][['segment_group']]
    # merge allowed segments:
    df38 = pd.concat([df37, df36s])
    
    # now we start again but calculate reads mapped per VGB/segment group 
    # and we will cut out the ones we don't want to report at the end 
    
    # so group by VGB: 
    df4 = df2.copy().groupby("VGB", as_index=False).sum()[['VGB','reads_mapped']]
    df5 = df2.copy().groupby(["VGB","cluster_member"], as_index=False).sum()[['VGB','len']].groupby("VGB", as_index=False).quantile(0.95)
    df6 = pd.merge(df4.copy(), df5, on='VGB', how='inner')
    df6['observed_RPK'] = df6['reads_mapped']/(df6['len']/1000)
    df6 = df6.copy()[['VGB','observed_RPK']]

    # now group again for segments:
    df7 = pd.merge(df6, ref.copy()[['VGB','segment_group','conserved_VGB']].drop_duplicates(), on='VGB')
    df8 = df7.copy()[~df7['segment_group'].str.contains("segment")][['segment_group','observed_RPK']].query("observed_RPK>0")
    # df8 is the non-segments final calculation
    
    df9 = df7.copy()[df7['segment_group'].str.contains("segment")]
    df91 = df9.query("conserved_VGB=='conserved'")
    df92 = df9.query("conserved_VGB=='not_conserved'").query("observed_RPK>0")
    df10 = pd.concat([df91,df92])
    df11 = df10[['observed_RPK','segment_group']].groupby("segment_group", as_index=False).mean().query("observed_RPK>0")

    # merge back together and drop VGBs not reaching 75% marker len:
    df12 = pd.concat([df8, df11])
    df13 = pd.merge(df12, df38, on='segment_group', how='inner')
    
    # format in HUMAnN like style:
    df14 = pd.merge(df13.copy(), taxref.copy()[['segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB']].drop_duplicates(), on='segment_group', how='left')
    df15 = df14.copy()
    df14['BAQLaVa VGB'] = df14['segment_group'] + "|nucleotide"
    df15['BAQLaVa VGB'] = df15['segment_group']
    df16 = pd.concat([df14, df15]).sort_values(by=['segment_group', 'BAQLaVa VGB']).rename(columns={'observed_RPK':newbase})[['BAQLaVa VGB',newbase, 'Taxonomy','Reference Species', 'Other ICTV Genomes in VGB']]

    # format tempfile to save out:
    df17 = df3.copy()[['# Gene Family','VGB', 'reads_mapped', 'marker_len_adj']]
    df17['RPK'] = df17['reads_mapped']/(df17['marker_len_adj']/1000)
    df18 = df17[df17['RPK']!=0][['# Gene Family','VGB','RPK']].rename(columns={"# Gene Family":"Marker"})
    return df16, df18

def process_baqlava_translated(base, format_df, ref, taxref):
    newbase = "_".join(base.split("_translated_"))
    df1 = format_df.copy()
    
    df2 = pd.merge(df1, ref.copy(), left_on='# Gene Family', right_on='protein', how='right')
    df2[base] = df2[base].fillna(0)
    
    #cut down to just the segment groups where at least one observation has been non-zero:
    df3 = df2.copy()[df2[base]!=0][['segment_group']].drop_duplicates()
    df4 = pd.merge(df2, df3, on='segment_group', how='inner')
    
    #turn df into dict to move through the averaging more quickly:
    dct = {}
    for i in range(len(df4)):
        if df4.iloc[i]['cluster_member'] in dct:
            dct[df4.iloc[i]['cluster_member']].append(df4.iloc[i][base])
        else:
            dct[df4.iloc[i]['cluster_member']] = [df4.iloc[i][base]]
    
    abund = []
    quantile = []
    for j in dct.keys():
        arr = np.array(dct[j])
        # 50:95 window:
        lo, hi = scipy.stats.mstats.mquantiles( arr, [.5, .95] )
        arr2 = arr[np.where((arr >= lo) & (arr <= hi))]
        abund.append(np.mean(arr2))
        # N% required to be detected:
        quant = scipy.stats.mstats.mquantiles( arr, [0.25] )
        quantile.append(quant[0])
        
    df5 = pd.DataFrame({'cluster_member':list(dct.keys()), 'observed_RPK':abund, 'quantile':quantile})
    df6 = pd.merge(df5, ref.copy()[['cluster_member','VGB','segment_group','conserved_VGB']].drop_duplicates(), on='cluster_member', how='left')

    # drop VGBs where <75% of proteins IDd within one cluster representative: 
    df61 = df6.copy()[df6['quantile']!=0]
    df62 = pd.merge(df6.copy().drop(columns={'observed_RPK'}), df61[['cluster_member','VGB', 'observed_RPK']], on=['cluster_member','VGB'], how='outer')
    df62['observed_RPK'] = df62['observed_RPK'].fillna(0)
    df63 = df62.sort_values(by='observed_RPK', ascending=False).groupby(["VGB"], as_index=False).first()
    
    # non-segmented:
    df7 = df63.copy()[~df63['segment_group'].str.contains("segment")].query("observed_RPK!=0").drop(columns=['conserved_VGB', 'VGB'])

    # segmented:
    df8 = df63.copy()[df63['segment_group'].str.contains("segment")]
    df91 = df8.copy().query("conserved_VGB=='conserved'")
    df92 = df8.copy().query("conserved_VGB=='not_conserved'").query("observed_RPK!=0")
    df9 = pd.concat([df91, df92])
    df10 = df9[['segment_group','observed_RPK']].groupby("segment_group", as_index=False).mean()

    # merge:
    df11 = pd.concat([df7, df10]).drop(columns={'quantile','cluster_member'}).query("observed_RPK!=0")
    
    # format in HUMAnN like style:
    df12 = pd.merge(df11.copy(), taxref.copy()[['segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB']].drop_duplicates(), on='segment_group', how='left')
    df13 = df12.copy()
    df13['BAQLaVa VGB'] = df13['segment_group'] + "|translated"
    df12['BAQLaVa VGB'] = df12['segment_group']
    df14 = pd.concat([df12, df13]).sort_values(by=['segment_group', 'BAQLaVa VGB']).rename(columns={'observed_RPK':newbase})[['BAQLaVa VGB',newbase, 'Taxonomy','Reference Species', 'Other ICTV Genomes in VGB']]
    
    # format tempfile to save out:
    df15 = df2.copy()[['# Gene Family','VGB', base]]
    df16 = df15[df15[base]!=0].rename(columns={"# Gene Family":"Protein", base:"RPK"})
    
    return df14, df16

def join_nuc_trans(nuc, trans, nucbase, ref):
    newbase = "_".join(nucbase.split("_nucleotide_"))
    
    nuc = nuc.rename(columns={newbase:'nucleotide'})
    nuc = nuc[~nuc['BAQLaVa VGB'].str.contains("\|")]
    trans = trans.rename(columns={newbase:'translated'})
    trans = trans[~trans['BAQLaVa VGB'].str.contains("\|")]
    df1 = pd.merge(nuc[['BAQLaVa VGB', 'nucleotide']], trans[['BAQLaVa VGB', 'translated']], on='BAQLaVa VGB', how='outer').fillna(0)
    
    # where we can, subtract abundance:
    df2 = df1.copy()[df1['translated']>df1['nucleotide']]
    df2 = df2.copy()
    df2['extra_protein'] = df2['translated']-df2['nucleotide']
    df3 = df1.copy()[~(df1['translated']>df1['nucleotide'])]
    
    df4 = pd.concat([df2, df3]).fillna(0)
    df4['Total'] = df4['nucleotide'] + df4['extra_protein']

    df5 = df4.melt(id_vars='BAQLaVa VGB').sort_values(by=['BAQLaVa VGB','variable'])
    df6 = df5.query("variable!='extra_protein'").query("value!=0")
    df7 = pd.merge(df6, ref.copy()[['segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB']].drop_duplicates(), left_on='BAQLaVa VGB', right_on='segment_group', how='left')
    
    # format similar to HUMAnN:
    df71 = df7.copy().query("variable=='Total'")
    df72 = df7.copy().query("variable!='Total'")
    df72['BAQLaVa VGB'] = df72['BAQLaVa VGB'] + "|" + df72['variable']

    if "bacterial_depleted" in nucbase:
        outcol = "_".join(nucbase.split("_bacterial_depleted_nucleotide_"))
        df8 = pd.concat([df71, df72]).rename(columns={'value':outcol}).reset_index().sort_values(by='index')[['BAQLaVa VGB',outcol,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']]
    else:
        df8 = pd.concat([df71, df72]).rename(columns={'value':newbase}).reset_index().sort_values(by='index')[['BAQLaVa VGB',newbase,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']]

    return df8


def run_reconciliation(sys1, sys2, sys3, nucref, transref, taxref, depl, inp_fa):
    #1: nuc only
    #2: trans only
    #3: both
    if sys1 == "1" or sys1 == "3":
        print('processing nuc')
        a,b = format_file_and_base(sys2)
        c = get_read_length(inp_fa)
        d,e = process_baqlava_nucleotide(a, b, nucref, c, taxref)
        e.to_csv(sys.argv[7], sep="\t", index=False)
    if sys1 == "2" or sys1 == "3":
        print('processing trans')
        f,g = format_file_and_base(sys3)
        h,i = process_baqlava_translated(f, g, transref, taxref)
        i.to_csv(sys.argv[8], sep="\t", index=False)
    if sys1 == "1":
        #if depl == "0":
        #    return d
        #elif depl == "1":
        #    d2 = d.rename(columns={a:a.replace("_bacterial_depleted","")})
        #    return d2
        return d
    elif sys1 == "2":
        #if depl == "0":
        #    return h
        #elif depl == "1":
        #    h2 = h.rename(columns={f:f.replace("_bacterial_depleted","")})
       	#    return h2
        return h
    elif sys1 == "3":
        j = join_nuc_trans(d, h, a, taxref)
        #if depl == "0":
        #    return j
        #elif depl == "1":
        #    j2 = j.rename(columns={a:a.replace("_bacterial_depleted","")})
        #    return j2
        return j
    else:
        return 'ERROR: incorrect reconciliation task requested'


def remove_problem_VGBs(df):
    # these are VGBs that contain bacterial genomes
    # this happened because ICTV reported accession as virus
    # that were bacterial genomes containing integrated virus
    df1 = df.copy()
    df1 = df1[~df1["BAQLaVa VGB"].str.contains("VGB_105795")]
    df1 = df1[~df1["BAQLaVa VGB"].str.contains("VGB_115203")]
    df1 = df1[~df1["BAQLaVa VGB"].str.contains("VGB_104822")]
    df1 = df1[~df1["BAQLaVa VGB"].str.contains("VGB_69885")]
    df1 = df1[~df1["BAQLaVa VGB"].str.contains("VGB_119070")]
    df1 = df1[~df1["BAQLaVa VGB"].str.contains("VGB_106952")]
    return df1

baq_out = run_reconciliation(sys.argv[1], sys.argv[2], sys.argv[3], nucleotide_reference, protein_reference, VGB_taxonomy, sys.argv[4], sys.argv[5])
baq_out2 = remove_problem_VGBs(baq_out)
baq_out2.to_csv(sys.argv[6], sep="\t", index=False)
