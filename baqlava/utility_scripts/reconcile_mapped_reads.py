import pandas as pd
import os
import sys
import numpy as np


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

marker_taxonomy = os.path.abspath(config.get('utility','marker_taxonomy'))
translated_marker_conversion = os.path.abspath(config.get('utility','translated_marker_conversion'))
species_conversion = os.path.abspath(config.get('utility','species_conversion'))

marker_taxonomy_df = pd.read_csv(marker_taxonomy, sep="\t", index_col='Unnamed: 0', low_memory=False)
translated_marker_conversion_df = pd.read_csv(translated_marker_conversion, sep="\t", index_col='Unnamed: 0')
species_conversion_df = pd.read_csv(species_conversion, sep="\t", index_col='Unnamed: 0')

def format_file_and_base(in_genefams):
    # get column of interest (basename plus RPK)
    base = ".".join(os.path.split( in_genefams )[1].split(".")[:-1]).replace('_genefamilies', '_Abundance-RPKs')
    # remove duplicate abundance lines (since everything is unclassified)
    df1 = pd.read_csv(in_genefams, sep="\t").copy()
    df2 = df1[~df1['# Gene Family'].str.contains("\|")]
    df3 = df2[df2['# Gene Family'] != 'UNMAPPED']
    return base, df3


def process_baqlava_nucleotide(base, format_df, nuc_tax, readlen):
    df2 = format_df.copy()
    df3 = df2.copy()[~df2['# Gene Family'].str.contains("_P1")]
    df3 = df3.copy()[~df3['# Gene Family'].str.contains("_P2")]
    df3 = df3.copy()[~df3['# Gene Family'].str.contains("_P3")]
    df3 = df3.copy()[~df3['# Gene Family'].str.contains("_M1")]
    df3 = df3.copy()[~df3['# Gene Family'].str.contains("_M2")]
    df3 = df3.copy()[~df3['# Gene Family'].str.contains("_M3")]
    
    df4 = pd.merge(df3.copy(), nuc_tax.copy(), right_on='marker', left_on='# Gene Family', how='inner').drop_duplicates()
    df4['reads_mapped'] = df4[base] * ((df4['marker_len']- (readlen-1))/1000) 
    df4['marker_len'] = df4['marker_len'].astype(float)
    df4 = df4[df4['marker_len']>400]
    df4['marker_len_adj'] = df4['marker_len'] - (readlen-1) 
    df5 = df4.copy().groupby("VGB", as_index=False).sum()[['VGB','reads_mapped']]
    df6 = df4.copy().groupby(["cluster_member","VGB"], as_index=False).sum()[['VGB','marker_len_adj']].groupby("VGB", as_index=False).quantile(0.95)
    df7 = pd.merge(df5, df6, on='VGB', how='inner')
    df7['observed_RPK'] = df7['reads_mapped']/(df7['marker_len_adj']/1000)
    df8 = df7.copy()[['VGB','observed_RPK']].drop_duplicates()
    
    #handle segmentation:
    df9 = nuc_tax.copy().dropna(subset=["segment_group"]).query("marker_len>400").query("marker!='no_marker'")[['Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species', 'VGB_specificity', 'segment_group', 'VGB_specificity_type', 'Virus name(s)', 'VGB']].drop_duplicates()
    df10 = pd.merge(df8, df9, on='VGB', how='inner')
    
    if len(df10) != 0:
        seg_grps = df10.copy()[['segment_group']].drop_duplicates()
        sdf1 = pd.merge(df9, seg_grps, on='segment_group', how='inner')
        sdf2 = pd.merge(sdf1, df10[['VGB', 'observed_RPK']], on='VGB', how='outer')
        sdf2['observed_RPK'] = sdf2['observed_RPK'].fillna(0)
        
        # drop specific VGB zeros and average remaining VGB abundances:
        sdf3 = sdf2.copy().query("VGB_specificity!='conserved'").query("VGB_specificity!='group1'").query("observed_RPK==0")
        sdf3['drop'] = 'yes'
        sdf4 = pd.merge(sdf2, sdf3[['VGB','drop']], on='VGB', how='outer').query("drop!='yes'")
        sdf5 = sdf4.copy().groupby("segment_group", as_index=False).mean()
        sdf6 = pd.merge(df9.copy()[['segment_group','Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']].drop_duplicates(), sdf5, on='segment_group', how='right')
        
        # merge back into VGB abundance table:
        df10['drop'] = 'yes'
        df11 = pd.merge(df10[['VGB','drop']], df8, on='VGB', how = 'outer').query("drop!='yes'")
        df12 = pd.merge(df11, nuc_tax.copy()[['VGB','Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']].drop_duplicates(), on='VGB', how='inner')
        df13 = pd.concat([df12, sdf6])
        df13['VGB'] = df13['VGB'].fillna(df13['segment_group'])
        return df13[['VGB', 'observed_RPK', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']]
    else:
        df12 = pd.merge(df8, nuc_tax.copy()[['VGB','Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']].drop_duplicates(), on='VGB', how='inner')
        return df12[['VGB', 'observed_RPK', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']]


def process_baqlava_translated(base, format_df, prot_conv, nuc_tax, readlen):
    df2 = format_df.copy()
    df31 = df2.copy()[df2['# Gene Family'].str.contains("_P1")]
    df32 = df2.copy()[df2['# Gene Family'].str.contains("_P2")]
    df33 = df2.copy()[df2['# Gene Family'].str.contains("_P3")]
    df34 = df2.copy()[df2['# Gene Family'].str.contains("_M1")]
    df35 = df2.copy()[df2['# Gene Family'].str.contains("_P2")]
    df36 = df2.copy()[df2['# Gene Family'].str.contains("_P3")]
    df4 = pd.concat([df31, df32, df33, df34, df35, df36])
    
    df6 = pd.merge(df4, prot_conv, left_on='# Gene Family', right_on='translated_marker', how='left')
    df7 = df6.copy().groupby("nucleotide_marker", as_index=False).sum()

    df8 = pd.merge(df7.copy(), nuc_tax.copy(), right_on='marker', left_on='nucleotide_marker', how='inner').drop_duplicates()
    df8['reads_mapped'] = df8[base] * ((df8['marker_len']- (readlen-1))/1000) # - n1 + 1
    df8['marker_len'] = df8['marker_len'].astype(float)
    df9 = df8[df8['marker_len']>400]
    df9 = df9.copy()
    df9['marker_len_adj'] = df9['marker_len'] - (readlen-1) 
    
    df10 = df9.copy().groupby("VGB", as_index=False).sum()[['VGB','reads_mapped']]
    df11 = df9.copy().groupby(["cluster_member","VGB"], as_index=False).sum()[['VGB','marker_len_adj']].groupby("VGB", as_index=False).quantile(0.95)
    df12 = pd.merge(df10, df11, on='VGB', how='inner')
    df12['observed_RPK'] = df12['reads_mapped']/(df12['marker_len_adj']/1000)
    df13 = df12.copy()[['VGB','observed_RPK']].drop_duplicates()
    
    #handle segmentation:
    df14 = nuc_tax.copy().dropna(subset=["segment_group"]).query("marker_len>400").query("marker!='no_marker'")[['Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species', 'VGB_specificity', 'segment_group', 'VGB_specificity_type', 'Virus name(s)', 'VGB']].drop_duplicates()
    df15 = pd.merge(df13, df14, on='VGB', how='inner')
    
    if len(df15) != 0:
        seg_grps = df15.copy()[['segment_group']].drop_duplicates()
        sdf1 = pd.merge(df14, seg_grps, on='segment_group', how='inner')
        sdf2 = pd.merge(sdf1, df15[['VGB', 'observed_RPK']], on='VGB', how='outer')
        sdf2['observed_RPK'] = sdf2['observed_RPK'].fillna(0)
        
        # drop specific VGB zeros and average remaining VGB abundances:
        sdf3 = sdf2.copy().query("VGB_specificity!='conserved'").query("VGB_specificity!='group1'").query("observed_RPK==0")
        sdf3['drop'] = 'yes'
        sdf4 = pd.merge(sdf2, sdf3[['VGB','drop']], on='VGB', how='outer').query("drop!='yes'")
        sdf5 = sdf4.copy().groupby("segment_group", as_index=False).mean()
        sdf6 = pd.merge(df14.copy()[['segment_group','Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']].drop_duplicates(), sdf5, on='segment_group', how='right')
        
        # merge back into VGB abundance table:
        df15['drop'] = 'yes'
        df16 = pd.merge(df15[['VGB','drop']], df13, on='VGB', how = 'outer').query("drop!='yes'")
        df17 = pd.merge(df16, nuc_tax.copy()[['VGB','Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']].drop_duplicates(), on='VGB', how='inner')
        df18 = pd.concat([df17, sdf6])
        df18['VGB'] = df18['VGB'].fillna(df18['segment_group'])
        df19 = df18[df18['observed_RPK']>=10]
        return df19[['VGB', 'observed_RPK', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']]
        
    else:
        df20 = pd.merge(df13, nuc_tax.copy()[['VGB','Realm','Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']].drop_duplicates(), on='VGB', how='inner')
        df21 = df20[df20['observed_RPK']>=10]
        return df21[['VGB', 'observed_RPK', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']] 

    
def merge_and_format_outputs(df1, df2, format_file, base):
    df1['Database'] = 'nucleotide'
    df1['sort'] = 2
    df2['Database'] = 'translated'
    df2['sort'] = 3
    df3 = pd.concat([df1[['VGB','Species','observed_RPK', 'Database','sort']], df2.rename(columns={'Genus_translated_collapse':'Genus'})[['VGB','Species','observed_RPK', 'Database', 'sort']]])
    df4 = pd.merge(format_file.copy(), df3, on=['VGB','Species'], how='right')
    df5 = df4[['Species_formatted','observed_RPK','Database','sort']].rename(columns={'Species_formatted':'Species'})
    df6 = df5.copy().groupby("Species", as_index=False).sum()
    df6['Database'] = "TOTAL"
    df6['sort'] = 1
    df7 = pd.concat([df5, df6]).sort_values(by = ['Species', 'sort'])
    df8 = df7.rename(columns={'observed_RPK':base})[['Species',base,'Database']]
    return df8
    
def run_functions_and_obtain_final_output(in_genefams, nuc_tax, prot_conv, readlen, format_file):
    a = format_file_and_base(in_genefams)
    b = process_baqlava_nucleotide(a[0], a[1], nuc_tax, readlen)
    c = process_baqlava_translated(a[0], a[1], prot_conv, nuc_tax, readlen)
    d = merge_and_format_outputs(b, c, format_file, a[0])
    return d


final_baq_profile = run_functions_and_obtain_final_output(sys.argv[1], marker_taxonomy_df, translated_marker_conversion_df, 100, species_conversion_df)


bas = os.path.split( sys.argv[1] )[1].split("_genefamilies.tsv")[0]
loc = os.path.split( sys.argv[1] )[0]
final_baq_profile.to_csv(loc + "/" + bas + "_BAQLaVa_profile.tsv", sep="\t", index=False)
