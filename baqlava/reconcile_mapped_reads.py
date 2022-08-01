import pandas as pd
import sys
import os

prot_tax_file = pd.read_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/additional_files/uniref90_ICTV_taxonomy_2.txt", sep="\t", index_col='Unnamed: 0', low_memory=False)
nuc_tax_file = pd.read_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/additional_files/nucleotide_ICTV_taxonomy_2.txt", sep="\t", index_col='Unnamed: 0', low_memory=False)

def process_baqlava_output(sysargv1, nuc_tax, prot_tax):
    #get column of interest (basename plus RPK)
    base = os.path.split( sysargv1 )[1].split(".")[0].replace('_genefamilies', '_Abundance-RPKs')
    #process nucleotide
    df1 = pd.read_csv(sysargv1, sep="\t").copy()
    df2 = df1[~df1['# Gene Family'].str.contains("\|")]
    df3 = pd.merge(df2.copy(), nuc_tax.copy(), right_on='genome', left_on='# Gene Family', how='inner')
    df3[['Species']] = df3['Species'].fillna(df3['sseqid'])
    df4 = df3.groupby(['Species'], as_index=False).sum()[['Species', base]]
    df4['subset'] = df4['Species'] + "|nucleotide"
    #process protein
    df5 = df2.copy()[df2['# Gene Family'].str.contains("UniRef90")]
    df6 = pd.merge(df5, prot_tax.copy(), right_on='UniRef90', left_on='# Gene Family', how='outer')
    df6[[base]] = df6[base].fillna(0)
    df7 = df6.groupby(['Species'], as_index=False).mean()[['Species', base]]
    df8 = df7[df7[base] != 0]
    df8['subset'] = df8['Species'] + "|translated"
    #now merge together and format results:
    df9 = pd.concat([df4,df8])
    df10 = df9.groupby(['Species'], as_index=False).sum()#[['Species', '1M_Abundance-RPKs']]
    df11 = pd.concat([df9, df10])
    df11[['subset']] = df11['subset'].fillna(df11['Species'])
    return df11.sort_values(by=['subset'], ascending=[True])[['subset', base]].rename(columns={'subset':'# Gene Family'})


processed_output = process_baqlava_output(sys.argv[1], nuc_tax_file, prot_tax_file)

bas = os.path.split( sys.argv[1] )[1].split("_genefamilies.tsv")[0]
loc = os.path.split( sys.argv[1] )[0]
processed_output.to_csv(loc + "/" + bas + "_baqlava_genefamilies.tsv", sep="\t", index=False)
