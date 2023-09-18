import pandas as pd
import sys

#sysargv1: final_databases/BAQLaVa_nucleotidedb_reference_file.txt
#sysargv2: BAQ_clu.tsv (mmseqs)
#sysargv3: bowtie2_clu.txt
#sysargv4: duplicated genomes
#sysargv5: connected components
#sysargv6: location to save out final df 

def combine_clustering(mm_clus_file, bt_clus_file, cc_file):
    mmclus = pd.read_csv(mm_clus_file, sep="\t", names=['cluster_rep_mm', 'cluster_member_mm'])
    btclus = pd.read_csv(bt_clus_file, sep="\t")
    # first join the clusters where the cluster rep from mmseqs2 (cluster 1) stayed a cluster rep in bowtie2 (cluster 2)
    clus = pd.merge(btclus, mmclus, left_on='cluster_member', right_on='cluster_rep_mm', how='inner')[['cluster_rep','cluster_member_mm']].drop_duplicates().rename(columns={'cluster_member_mm':'cluster_member'})
    clus['clus_source'] = 'bowtie2'
    # code not needed since all members are included in first merge
    #clus = pd.concat([btclus, additional_members]).drop_duplicates()
    # add back in the mmseqs2 genomes that were not clustered further by bowtie2:
    clus2 = pd.merge(btclus, mmclus, left_on='cluster_member', right_on='cluster_rep_mm', how='outer', indicator = True).query("_merge=='right_only'")[['cluster_member_mm','cluster_rep_mm']].rename(columns={'cluster_member_mm':'cluster_member','cluster_rep_mm':'cluster_rep'})
    clus2['clus_source'] = 'mmseqs2'
    clus3 = pd.concat([clus,clus2]).drop_duplicates()
    # check that the new clustered df represents the correct number of genomes:
    if len(clus3) == len(mmclus):
        ccs = pd.read_csv(cc_file, sep='\t', index_col='Unnamed: 0').query("genomes!='*'")
        clus4 = pd.merge(ccs, clus3, left_on='genomes', right_on='cluster_member', how='inner')
        if len(clus4) == len(clus3):
            clus4['connected_component'] = clus4['connected_component'].astype(str)
            clus4['VGB'] = 'VGB_' + clus4['connected_component']
            return clus4[['VGB','connected_component','cluster_rep','cluster_member','clus_source']]
        else: return "ERROR"
    else:
        return "ERROR"
        
def remove_dups_for_final_clusters(clus_df, ref_file, dup_file):
    ref = pd.read_csv(ref_file, sep="\t", index_col='Unnamed: 0')[['Orig_Header','New_Header','database','length']]
    df1 = pd.merge(clus_df, ref, left_on='cluster_member', right_on='New_Header', how='inner')
    dups = pd.read_csv(dup_file, sep="\t", index_col='Unnamed: 0').rename(columns={'genome':'cluster_member'})
    dups['duplicated'] = 'yes'
    df2 = pd.merge(df1, dups, on='cluster_member', how='outer')
    df2['duplicated'] = df2['duplicated'].fillna('no')
    
    # get the connected components where cluster reps were duplicated genomes:
    # we would like to throw these out, but only when doing so doesn't remove all ability to markerize 
    # (eg if all cluster reps are duplicates - then we would throw out all 100mers and can't markerize)
    ccs_with_duplicated_clusterreps =  set(df2.groupby(['connected_component','cluster_rep','duplicated'], as_index=False).first().query("duplicated=='yes'")['connected_component'])
    ccs_with_nondup__clusterreps =  set(df2.groupby(['connected_component','cluster_rep','duplicated'], as_index=False).first().query("duplicated=='no'")['connected_component'])
    ccs_with_duplicated_clusterreps_to_remove = ccs_with_duplicated_clusterreps.intersection(ccs_with_nondup__clusterreps)
    # merge back in the info for ccs where we want to drop duplicates
    df3 = pd.DataFrame({'connected_component':list(ccs_with_duplicated_clusterreps_to_remove)})
    df3['drop_duplicated'] = 'yes'
    df4 = pd.merge(df3, df2, on='connected_component', how='outer')
    df4['drop_duplicated'] = df4['drop_duplicated'].fillna('no')
    # remove duplciates:
    # keep only lines that don't meet the condition drop_duplicated=='yes' AND duplicated=='yes'
    df5 = df4.query("drop_duplicated=='no'")
    df6 = df4.query("drop_duplicated=='yes'").query("duplicated=='no'")
    df7 = pd.concat([df5, df6])
    df7['connected_component'] = df7['connected_component'].astype(str)
    df7['VGB'] = 'VGB_' + df7['connected_component']

    return df7        

cluster1 = combine_clustering(sys.argv[2], sys.argv[3], sys.argv[5]).to_csv(sys.argv[6]+"final_clusters.txt", sep="\t")
#remove_dups_for_final_clusters(cluster1, sys.argv[1], sys.argv[4]).to_csv(sys.argv[6]+"final_clusters.txt", sep="\t")
