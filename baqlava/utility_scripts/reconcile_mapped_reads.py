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

nucleotide_taxonomy = os.path.abspath(config.get('utility','nucleotide_taxonomy'))
protein_taxonomy = os.path.abspath(config.get('utility','protein_taxonomy'))
viral_name_conversion = os.path.abspath(config.get('utility','viral_name_conversion'))
jaccard = os.path.abspath(config.get('utility','jaccard'))
segmented = os.path.abspath(config.get('utility','segmented'))

nuc_tax_file = pd.read_csv(nucleotide_taxonomy, sep="\t", index_col='Unnamed: 0', low_memory=False)
prot_tax_file = pd.read_csv(protein_taxonomy, sep="\t", index_col='Unnamed: 0', low_memory=False)
                           
viral_names2 = pd.read_csv(viral_name_conversion, sep="\t", index_col='Unnamed: 0', low_memory=False)

jaccard_df = pd.read_csv(jaccard, sep="\t", low_memory=False)
segment_ref = pd.read_csv(segmented, sep="\t", index_col='Unnamed: 0')

def format_file_and_base(sysargv1):
    #get column of interest (basename plus RPK)
    base = ".".join(os.path.split( sysargv1 )[1].split(".")[:-1]).replace('_genefamilies', '_Abundance-RPKs')
    #process nucleotide
    df1 = pd.read_csv(sysargv1, sep="\t").copy()
    df2 = df1[~df1['# Gene Family'].str.contains("\|")]
    df3 = df2[df2['# Gene Family'] != 'UNMAPPED']
    return base, df3

def process_baqlava_nucleotide_step1(base, format_df, nuc_tax, segment_df):
    df2 = format_df.copy()
    df3 = df2.copy()[~df2['# Gene Family'].str.contains("UniRef90")]
    #process nucleotide hits
    df4 = pd.merge(df3.copy(), nuc_tax.copy(), right_on='sseqid', left_on='# Gene Family', how='inner').drop_duplicates()
    df5 = pd.merge(df4, segment_df, on=['Genus','Species'], how='outer', indicator = True)
    df6 = df5.query("_merge=='both'") # segemented genomes
    df7 = df5.query("_merge=='left_only'") # all other hits to unknowns or non segemented genomes
    return df6, df7
    
def process_baqlava_nucleotide_step2_segmented(base, step1_df, segment_df):    
    df1 = step1_df.copy()
    #filter genomes w/o 50% segments present
    df2 = df1.groupby(['Genus','Species'], as_index=False).count()[['Genus','Species','# segments']].rename(columns={'# segments':'# segments observed'})
    df3 = pd.merge(df2, segment_df, on=['Genus','Species'], how='inner')
    df3['% segments observed'] = df3['# segments observed'] / df3['# segments']
    df4 = df3[df3['% segments observed'] >= 0.5]
    df5 = pd.merge(df4.copy()[['Genus','Species']], df1.copy(), on=['Genus','Species'], how='inner').groupby(['Genus','Species'], as_index=False).sum()[['Genus','Species',base]]

    df6 = pd.merge(df4.copy()[['Genus','Species','# segments']], df5, on=['Genus','Species'], how='inner')
    df6['segment_adjusted_abundance'] = df6[base] / df6['# segments']
    #great! now format it all back nicely so we can merge with non-segmented hits
    df7 = df6[['Genus', 'Species', 'segment_adjusted_abundance']].rename(columns={'segment_adjusted_abundance':base})
    return df7 

def process_baqlava_nucleotide_step3_nonsegmented(base, step1_df):  
    df1 = step1_df.copy()
    df1[['Species']] = df1['Species'].fillna(df1['sseqid'])
    df1['Genus'] = df1['Genus'].fillna('unknown')
    df2 = df1.groupby(['Genus', 'Species'], as_index=False).sum()[['Genus', 'Species', base]]
    return df2

def process_baqlava_nucleotide_step4_combine_and_format(base, seg_df, nonseg_df, names_file): 
    df1 = pd.concat([seg_df.copy(), nonseg_df.copy()])
    df2 = pd.merge(df1, names_file, on='Species', how='left').drop_duplicates()    
    return df2

def process_baqlava_translated_step1(base, format_df, prot_tax):
    df2 = format_df.copy()
    df3 = df2.copy()[df2['# Gene Family'].str.contains("UniRef90")]
    df4 = pd.merge(df3, prot_tax.copy(), right_on='UniRef90', left_on='# Gene Family', how='outer')
    df4[[base]] = df4[base].fillna(0)
    return df4

def process_baqlava_translated_step2(base, step1_df):
    df1 = step1_df.copy()
    df2 = df1[['UniRef90','Species','# Gene Family']]
    df3 = df2.groupby('Species', as_index=False).count().rename(columns={'UniRef90':'#_UR90s_species','# Gene Family':'#_UR90s_observed'})
    df3['percent_UR90s_observed'] = df3['#_UR90s_observed']/df3['#_UR90s_species']
    df4 = df3[df3['percent_UR90s_observed']>= 0.5]
    df5 = pd.merge(df1, df4[['Species']], on='Species', how='inner')
    return df4, df5
    
def process_baqlava_translated_step3(base, step2_df, percent_observed_df):  
    retdf = pd.DataFrame({'Species':[], 'UniRef90':[], base:[]})
    df1 = step2_df.copy()
    specs = list(set(df1['Species']))
    for i in specs:
        tempdf = df1[df1['Species']==i][['Species','UniRef90',base]]
        hi_window = int(len(tempdf) * 0.95)
        lo_window = int(len(tempdf) * 0.75)
        tempdf = tempdf.sort_values(by=[base])
        tempdf = tempdf.iloc[lo_window:hi_window+1]
        retdf = pd.concat([retdf,tempdf])
    df2 = retdf.groupby('Species', as_index=False).mean()
    df3 = pd.merge(df2, percent_observed_df.copy(), on='Species', how='inner')
    return df3[['Species',base,'percent_UR90s_observed']].drop_duplicates()

def process_baqlava_translated_step4(base, step3_df, jacc_df):
    df1 = step3_df.copy()
    df2 = df1.sort_values(by=['percent_UR90s_observed',base], ascending = [False,False])
    skip = []
    #gen_list = []
    spec_list = []
    RPK_list = []
    clus_list = []
    for i in df2['Species']:
        if i in skip:
            pass
        else:
            if df2[df2['Species']==i].iloc[0][base] == 0:
                pass
            else:
                clus1 = jacc_df[jacc_df['Pair1']==i]
                clus2 = jacc_df[jacc_df['Pair2']==i]
                clus3 = pd.concat([clus1,clus2])
                clus3 = clus3[clus3['Jaccard']>0.05]
                if len(clus3) == 0:
                    spec_list.append(i)
                    RPK_list.append(df2[df2['Species']==i].iloc[0][base])
                    clus_list.append('no')
                else:
                    clus4 = pd.DataFrame({'Species':list(set(clus3['Pair1'])) + list(set(clus3['Pair2']))}).drop_duplicates()
                    clus5 = pd.merge(clus4, df2.copy(), on='Species', how='left').drop_duplicates().sort_values(by=['percent_UR90s_observed',base], ascending = [False,False])
                    hit = clus5.iloc[0]['Species']
                    spec_list.append(clus5.iloc[0]['Species'])
                    RPK_list.append(clus5.iloc[0][base])
                    clus_list.append('yes')
                    for j in clus5['Species']:
                        skip.append(j) 
    return pd.DataFrame({'Species':spec_list, base:RPK_list, 'cluster?':clus_list})

def format_nuc_trans_hits_for_return(base, nuc_abund_df, trans_abund_df, prot_tax):
    #format genus hits :
    prot_genus_names = prot_tax.copy()[['Genus','Species']].drop_duplicates().fillna('unknown')
    df1 = trans_abund_df.copy()
    df2 = pd.merge(df1, prot_genus_names, on='Species',how='left')
    df3 = df2.groupby('Genus', as_index=False).sum()
    df3['Virus'] = df3['Genus'] + "|unclassified"
    df3['Database'] = 'translated'
    df3['Virus_metadata'] = np.nan
    #format nuc hits :
    df4 = nuc_abund_df.copy()
    df4['Virus'] = df4['Genus'] + "|" +  df4['Species']
    df4 = df4.rename(columns={'Virus name(s)':'Virus_metadata'})
    df4['Database'] = 'nucleotide'
    return pd.concat([df3,df4[['Genus', base, 'Virus', 'Database', 'Virus_metadata']]]).sort_values(by=['Genus', 'Database'])[['Virus',base,'Virus_metadata','Database']]

def format_nuc_only_hits_for_return(base, nuc_abund_df):
    #format nuc hits :
    df4 = nuc_abund_df.copy()
    df4['Virus'] = df4['Genus'] + "|" +  df4['Species']
    df4 = df4.rename(columns={'Virus name(s)':'Virus_metadata'})
    df4['Database'] = 'nucleotide'
    return df4[['Genus', base, 'Virus', 'Database', 'Virus_metadata']].sort_values(by=['Genus', 'Database'])[['Virus',base,'Virus_metadata','Database']]

def run_functions_and_obtain_final_output(sys_argv_1, nuc_tax, prot_tax, segment_df, names_file, jacc_df):
    a = format_file_and_base(sys_argv_1)
    b, c = process_baqlava_nucleotide_step1(a[0], a[1], nuc_tax, segment_df)
    d = process_baqlava_nucleotide_step2_segmented(a[0], b, segment_df) 
    e = process_baqlava_nucleotide_step3_nonsegmented(a[0], c)
    f = process_baqlava_nucleotide_step4_combine_and_format(a[0], d, e, names_file)
    g = process_baqlava_translated_step1(a[0], a[1], prot_tax) #all protein sets per species
    h, i = process_baqlava_translated_step2(a[0], g) 
    if len(h) == 0:
        print("nuc only")
        j = format_nuc_only_hits_for_return(a[0], f).drop_duplicates()
        return j
    else:
        k = process_baqlava_translated_step3(a[0], i, h) #calculate 75:95 window abundance
        l = process_baqlava_translated_step4(a[0], k, jacc_df) #sort potential hits by jaccard similarity, remove FP
        m = format_nuc_trans_hits_for_return(a[0], f, l, prot_tax).drop_duplicates()
        return m

final_baq_profile = run_functions_and_obtain_final_output(sys.argv[1], nuc_tax_file, prot_tax_file, segment_ref, viral_names2, jaccard_df)

bas = os.path.split( sys.argv[1] )[1].split("_genefamilies.tsv")[0]
loc = os.path.split( sys.argv[1] )[0]
final_baq_profile.to_csv(loc + "/" + bas + "_BAQLaVa_profile.tsv", sep="\t", index=False)
