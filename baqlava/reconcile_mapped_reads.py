import pandas as pd
import sys

#load reference files:
nuc_ref_file = pd.read_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/baclava_nucleotide_annotations.txt", sep="\t", index_col='Unnamed: 0').rename(columns={'contig_id':'sseqid', 'name':'genome','length':'subject_length'})

prot_tax_file = pd.read_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/uniref90_ICTV_taxonomy.txt", sep="\t", index_col='Unnamed: 0')
nuc_tax_file = pd.read_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/nucleotide_ICTV_taxonomy.txt", sep="\t", index_col='Unnamed: 0')

#load in data:
HUMAnN_bt2 = pd.read_csv(sys.argv[1],
                         sep="\t",
                         names=['reads', 'hit','pident','alignment length','unknow3','unknow4','unknow5','unknow6','unknow7','unknow8','unknow9','unknow10'])[[
    'reads', 'hit','pident','alignment length','unknow6','unknow7','unknow8']]



HUMAnN_diam = pd.read_csv(sys.argv[2],
                         sep = "\t",
                          names = ['qseqid', 'sseqid','pident','alignment length','mismatch','gap open','qstart','qend','sstart','send','evalue','bitscore'])


def calculate_RPK_cov(bt2_df, ref_df, taxonomy):
    df1 = bt2_df.copy()
    df2 = df1.groupby(['hit'], as_index=False).count()
    ### inner below should be changes to 'left' in final code where ref file matches database mapped to I think
    df3 = pd.merge(df2, ref_df.copy(), left_on='hit', right_on='sseqid', how='inner')
    df3['RPK'] = (df3['reads']*1000)/df3['subject_length']
    df4 = pd.merge(df3, taxonomy.copy(), on='sseqid', how='inner')
    return df4[['hit', 'RPK', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum',
       'Class', 'Subclass', 'Order', 'Suborder', 'Family', 'Subfamily',
       'Genus', 'Subgenus', 'Species', 'genome']]


def format_diam_df(diam_df):
    df1 = diam_df.copy()
    df1[['qseqid', 'qlen']] = df1['qseqid'].str.split("|", expand=True)
    df1['qlen'] = df1.qlen.astype(int)
    df1['qseqid'] = df1.qseqid.astype(str)
    df1['sseqid'] = df1.sseqid.astype(str)
    print("part 1 done")
    df1['qcov'] = (abs(df1['qend'] - df1['qstart'])+1)/df1['qlen']
    df1['waafle'] = df1['qcov'] * df1['pident']
    print("part 2 done")
    df2 = df1.sort_values(by=['qseqid', 'waafle'], ascending=[True, False])
    df3 = df2.groupby(['qseqid'], as_index=False).first()
    df4 = df3.groupby(['sseqid'], as_index=False).count()[['sseqid', 'qseqid']].rename(columns={'qseqid':'mapped_reads'})
    return df4

def ur90_to_tax_RPK(diam_df):
    df1 = diam_df.copy()
    df1[['sseqid', 'slen']] = df1['sseqid'].str.split("|", expand=True)
    df1['slen'] = df1.slen.astype(int)
    df1['RPK'] = (df1['mapped_reads']*1000)/df1['slen']
    return df1

def format_prot_ref(df):
    df1 = df.copy()
    df1[['UniRef90', 'Length']] = df1['UniRef90'].str.split("|", expand=True)
    df1['Length'] = df1.Length.astype(int)
    return df1

def ur90_to_tax(diam_df, ur90_ref_df):
    df1 = pd.merge(diam_df.copy(), ur90_ref_df.copy(), left_on=['sseqid', 'slen'], right_on=['UniRef90', 'Length'], how='inner')
    df2 = df1.groupby(['Species'], as_index=False).mean()
    corr_tax = ur90_ref_df[['Realm', 'Subrealm', 'Kingdom',
       'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order',
       'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species']].drop_duplicates()
    df3 = pd.merge(df2[['Species', 'mapped_reads', 'Length', 'RPK']], corr_tax, on='Species', how='left')
    return df3[df3['RPK']>10]

#format protein taxonomy file for use:
prot_tax_file = format_prot_ref(prot_tax_file)

#process nucleotide mapped reads:
bowtie_mapping = calculate_RPK_cov(HUMAnN_bt2, nuc_ref_file, nuc_tax_file)

#process protein mapped reads:
diamond_mapping1 = format_diam_df(HUMAnN_diam)
diamond_mapping2 = ur90_to_tax_RPK(diamond_mapping1)
diamond_mapping3 = ur90_to_tax(diamond_mapping2, prot_tax_file)

#export data:
bowtie_mapping.to_csv("baqlava_nucleotide_profile.txt")
diamond_mapping3.to_csv("baqlava_protein_profile.txt")
