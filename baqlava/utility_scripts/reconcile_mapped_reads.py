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

nucleotide_reference = pd.read_csv(nucleotide_reference_file, sep="\t", index_col='Unnamed: 0', low_memory=False)
protein_reference = pd.read_csv(protein_reference_file, sep="\t", index_col='Unnamed: 0', low_memory=False)



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


def process_baqlava_nucleotide(base, format_df, ref, readlen):
    df1 = format_df.copy()
    
    df2 = pd.merge(df1, ref.copy(), left_on='# Gene Family', right_on='marker', how='right')
    df2[base] = df2[base].fillna(0)    
    df2['reads_mapped'] = df2[base] * ((df2['len']- (readlen-1))/1000) 
    df3 = df2[df2['len']>400]
    df3 = df3.copy()
    df3['marker_len_adj'] = df3['len'] - (readlen-1) 
    df4 = df3.copy().groupby("segment_group", as_index=False).sum()[['segment_group','reads_mapped']]
    df5 = df3.copy().groupby(["segment_group","cluster_member"], as_index=False).sum()[['segment_group','marker_len_adj']].groupby("segment_group", as_index=False).quantile(0.95)
    df6 = pd.merge(df4, df5, on='segment_group', how='inner')
    df6['observed_RPK'] = df6['reads_mapped']/(df6['marker_len_adj']/1000)
    df7 = df6.copy()[['segment_group','observed_RPK']].drop_duplicates().query("observed_RPK>0")
    
    df8 = pd.merge(df7, ref.copy()[['segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB']].drop_duplicates(), on='segment_group', how='left')
    
    # format in HUMAnN like style:
    df9 = df8.copy()
    df9['BAQLaVa VGB'] = df9['segment_group'] + "|nucleotide"
    df8['BAQLaVa VGB'] = df8['segment_group']
    df10 = pd.concat([df9, df8]).sort_values(by=['segment_group', 'BAQLaVa VGB']).rename(columns={'observed_RPK':base})

    # format tempfile to save out:
    df11 = df3.copy()[['# Gene Family','VGB', 'reads_mapped', 'marker_len_adj']]
    df11['RPK'] = df11['reads_mapped']/(df11['marker_len_adj']/1000)
    df12 = df11[df11['RPK']!=0][['# Gene Family','VGB','RPK']].rename(columns={"# Gene Family":"Marker"})
    return df10, df12


def process_baqlava_translated(base, format_df, ref):
    df1 = format_df.copy()
    
    df2 = pd.merge(df1, ref.copy(), left_on='# Gene Family', right_on='protein', how='right')
    df2[base] = df2[base].fillna(0)  
    
    #cut down to just the segment groups where at least one observation has been non-zero:
    df3 = df2.copy()[df2[base]!=0][['segment_group']].drop_duplicates()
    df4 = pd.merge(df2, df3, on='segment_group', how='inner')
    
    #turn df into dict to move through the averaging more quickly:
    dct = {}
    for i in range(len(df4)):
        if df4.iloc[i]['VGB'] in dct:
            dct[df4.iloc[i]['VGB']].append(df4.iloc[i][base])
        else:
            dct[df4.iloc[i]['VGB']] = [df4.iloc[i][base]]
    
    abund = []
    for j in dct.keys():
        arr = np.array(dct[j])
        lo, hi = scipy.stats.mstats.mquantiles( arr, [.5, .95] )
        arr2 = arr[np.where((arr >= lo) & (arr <= hi))]
        abund.append(np.mean(arr2))
    df5 = pd.DataFrame({'VGB':list(dct.keys()), 'observed_RPK':abund})
    df6 = pd.merge(df5, ref.copy()[['VGB','segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB', 'conserved_VGB']].drop_duplicates(), on='VGB', how='left')
    # non-segmented:
    df7 = df6.copy()[~df6['segment_group'].str.contains("segment")].query("observed_RPK!=0").drop(columns=['conserved_VGB', 'VGB'])
    
    # segmented:
    df8 = df6.copy()[df6['segment_group'].str.contains("segment")]
    df91 = df8.copy().query("conserved_VGB=='conserved'")
    df92 = df8.copy().query("conserved_VGB=='not_conserved'").query("observed_RPK!=0")
    df9 = pd.concat([df91, df92])
    df10 = df9[['segment_group','observed_RPK']].groupby("segment_group", as_index=False).mean()
    df11 = pd.merge(df10, ref.copy()[['segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB']].drop_duplicates(), on='segment_group', how='left')
    
    df12 = pd.concat([df7, df11])
    
    # format in HUMAnN like style:
    df13 = df12.copy()
    df13['BAQLaVa VGB'] = df13['segment_group'] + "|translated"
    df12['BAQLaVa VGB'] = df12['segment_group']
    df14 = pd.concat([df12, df13]).sort_values(by=['segment_group', 'BAQLaVa VGB']).rename(columns={'observed_RPK':base})

    # format tempfile to save out:
    df15 = df2.copy()[['# Gene Family','VGB', base]]
    df16 = df15[df15[base]!=0].rename(columns={"# Gene Family":"Protein", base:"RPK"})
    
    return df14, df16


def join_nuc_trans(nuc, trans, nucbase, transbase, ref):
    nuc = nuc.rename(columns={nucbase:'nucleotide'}).drop(columns=['BAQLaVa VGB']).drop_duplicates()
    trans = trans.rename(columns={transbase:'translated'}).drop(columns=['BAQLaVa VGB']).drop_duplicates()
    df1 = pd.merge(nuc[['segment_group', 'nucleotide']], trans[['segment_group', 'translated']], on='segment_group', how='outer').fillna(0)
    
    # where we can, subtract abundance:
    df2 = df1[df1['translated']>df1['nucleotide']]
    df2 = df2.copy()
    df2['extra_protein'] = df2['translated']-df2['nucleotide']
    df3 = df1[~(df1['translated']>df1['nucleotide'])]
    df4 = pd.concat([df2, df3]).fillna(0)
    df4['Total'] = df4['nucleotide'] + df4['extra_protein']
    df5 = df4.melt(id_vars='segment_group').sort_values(by=['segment_group','variable'])
    df6 = df5.query("variable!='extra_protein'").query("value!=0")
    df7 = pd.merge(df6, ref.copy()[['segment_group','Taxonomy','Reference Species','Other ICTV Genomes in VGB']].drop_duplicates(), on='segment_group', how='left')
    
    # format similar to HUMAnN:
    df71 = df7.copy().query("variable=='Total'")
    df71['BAQLaVa VGB'] = df71['segment_group']
    df72 = df7.copy().query("variable!='Total'")
    df72['BAQLaVa VGB'] = df72['segment_group'] + "|" + df72['variable']
    
    df8 = pd.concat([df71, df72]).rename(columns={'value':nucbase}).reset_index().sort_values(by='index')
    return df8[['BAQLaVa VGB',nucbase,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']]


def run_reconciliation(sys1, sys2, sys3, nucref, transref, depl, inp_fa):
    #1: nuc only
    #2: trans only
    #3: both
    if sys1 == "1" or sys1 == "3":
        print('processing nuc')
        a,b = format_file_and_base(sys2)
        c = get_read_length(inp_fa)
        d,e = process_baqlava_nucleotide(a, b, nucref, c)
        e.to_csv(sys.argv[7], sep="\t", index=False)
    if sys1 == "2" or sys1 == "3":
        print('processing trans')
        f,g = format_file_and_base(sys3)
        h,i = process_baqlava_translated(f, g, transref)
        i.to_csv(sys.argv[8], sep="\t", index=False)
    if sys1 == "1":
        if depl == "0":
            return d[['BAQLaVa VGB',a,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']]
        elif depl == "1":
            d2 = d[['BAQLaVa VGB',a,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']].rename(columns={a:a.replace(".bacterial_depleted_nucleotide","")})
            return d2
    elif sys1 == "2":
        if depl == "0":
            return h[['BAQLaVa VGB',f,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']]
        elif depl == "1":
            h2 = h[['BAQLaVa VGB',f,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']].rename(columns={f:f.replace(".bacterial_depleted_translated","")})
       	    return h2
    elif sys1 == "3":
        j = join_nuc_trans(d, h, a, f, nucref)
        if depl == "0":
            return j[['BAQLaVa VGB',a,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']]
        elif depl == "1":
            j2 = j[['BAQLaVa VGB',a,'Reference Species','Taxonomy','Other ICTV Genomes in VGB']].rename(columns={a:a.replace(".bacterial_depleted_nucleotide","")})
            return j2
    else:
        return 'ERROR: incorrect reconciliation task requested'
    
baq_out = run_reconciliation(sys.argv[1], sys.argv[2], sys.argv[3], nucleotide_reference, protein_reference, sys.argv[4], sys.argv[5])
baq_out.to_csv(sys.argv[6], sep="\t", index=False)
