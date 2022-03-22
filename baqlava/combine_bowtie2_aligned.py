import pandas as pd
import sys

HUMAnN_bt2 = pd.read_csv(sys.argv[1], sep="\t", names=['reads', 'hit','pident','alignment length','unknow3','unknow4','unknow5','unknow6','unknow7','unknow8','unknow9','unknow10'])
Viral_bt2 = pd.read_csv(sys.argv[2], sep="\t", names=['reads', 'hit','pident','alignment length','unknow3','unknow4','unknow5','unknow6','unknow7','unknow8','unknow9','unknow10'])

def get_saveout_name(sysargv1):
    name = str(sysargv1)
    name = name.split("/")[-1]
    name = name.replace('_bowtie2_aligned.tsv', '')
    name = name + '_combined_bowtie2_aligned.tsv'
    return name

def pull_reads_from_viral_mapping(df_H, df_V):
    HUMAnN_reads = set(df_H['reads'])
    viral_reads = set(df_V['reads'])
    unique_viral_reads = viral_reads - HUMAnN_reads
    unique_df = pd.DataFrame({'reads':list(unique_viral_reads)})
    viral_culled = pd.merge(unique_df, df_V, on='reads', how='left')
    combined_df = pd.concat([df_H,viral_culled])
    if len(combined_df) == len(viral_culled) + len(df_H):
        return combined_df
    else:
        print("error")

name = get_saveout_name(sys.argv[1])

pull_reads_from_viral_mapping(HUMAnN_bt2, Viral_bt2).to_csv(name, sep="\t", index=False, header=False)
