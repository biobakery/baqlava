import pandas as pd
import sys
import networkx as nx

#sysargv1 = mapped_reads file, only columns 1,2,3
#sysargv2 = "test/bowtie2/BAQLaVa_dereplicated_100mers.fasta"
#sysargv3 = "test/final_databases/BAQLaVa_nucleotidedb_reference_file.txt"
#sysargv4 = mmseqs clusters
#sysargv5 = saving location


def collect_all_possibe_markers(file):
    counter = 0
    genomes = []
    reads = []
    with open(file, "r") as fasta:
        for i in fasta:
            i = i.strip()
            if i[0] == ">":
                counter += 1
                genomes.append(i.split("_")[0][1:])
                reads.append(i.split("_")[1])
    return pd.DataFrame({'genome':genomes, 'read':reads}).groupby('genome', as_index=False).count()

def calculate_cov(markerdf, mappedreads_file):
    mapped1 = pd.read_csv(mappedreads_file, sep="\t")
    mapped2 = pd.merge(mapped1, markerdf, left_on='query', right_on='genome', how='left')
    mapped2['q_cov'] = mapped2['n_markers_mapped']/mapped2['read']
    mapped2 = mapped2.rename(columns={'read':'q_reads'})[['query','subject','n_markers_mapped','q_reads','q_cov']]
    mapped3 = mapped2[mapped2['q_cov']>=0.50]
    return mapped3

def build_connected_components(cov_df, mmseqs_file):
    
    df1 = pd.read_csv(mmseqs_file, sep="\t", names=['cluster_member','cluster_rep'])
    
    G = nx.Graph()
    G.add_nodes_from(list(df1['cluster_rep']))
    
    for i in range(len(df1)):
        G.add_edge(df1.iloc[i]['cluster_member'], df1.iloc[i]['cluster_rep'])
    
    for j in range(len(cov_df)):
        G.add_edge(cov_df.iloc[j]['query'], cov_df.iloc[j]['subject'])
    edges_list = []
    ccs = []
    counter = 0
    for cc in nx.connected_components(G):
        counter += 1
        for node in cc:
            edges_list.append(node)
            ccs.append(counter)
    return pd.DataFrame({'connected_component':ccs, 'genomes':edges_list})   



def cluster(markerdf, df, ref, loc):

    df1 = df.copy()
    df1 = df1[df1['subject'] != '*']
    df1 = df1[df1['q_cov']>0.85]
    df2 = pd.merge(df1, markerdf[['genome']], left_on='subject', right_on='genome', how='inner')[['query', 'subject', 'n_markers_mapped','q_reads', 'q_cov']]
    
    ref1 = pd.read_csv(ref, sep="\t", index_col='Unnamed: 0')[['New_Header','length','database']]
    df3 = pd.merge(df2, ref1, left_on='query', right_on='New_Header', how='left').rename(columns={'length':'q_length','database':'q_database'})[['query','subject','q_cov','q_length','q_database']]
    df4 = pd.merge(df3, ref1, left_on='subject', right_on='New_Header', how='left').rename(columns={'length':'s_length','database':'s_database'})[['query','subject','q_cov','q_length', 's_length','q_database','s_database']]
    df4['cluster_assigned'] = 'bowtie2'

    genomes_bylength = list(df4[['subject','s_length']].drop_duplicates().sort_values(by='s_length', ascending=False)['subject'])
    working_df = df4.copy()
    clustered_df = pd.DataFrame({'cluster_member':[],'cluster_rep':[]})
    clustered = set()
    counter = 0
    remove = []
    
    with open(loc + "bowtie2_clu.txt", "w") as clu:
        print("file open")
        clu.write("cluster_member")
        clu.write("\t")
        clu.write("cluster_rep")
        clu.write("\n")
        print("header written")
        for i in genomes_bylength:
            counter += 1
            if counter % 1000 == 0:
                print(counter)
            if i in clustered:
                pass
            else:
                tempdf = working_df.copy()[working_df['subject']==i]
        
                if len(tempdf) == 0:
                    pass
                else:
                    clu.write(i)
                    clu.write("\t")
                    clu.write(i)
                    clu.write("\n")
                    for j in tempdf['query']:
                        clu.write(j)
                        clu.write("\t")
                        clu.write(i)
                        clu.write("\n")
            
                    remove = list(tempdf['query'])
                    remove.append(i)

                    for j in remove:
                        working_df = working_df[working_df['query'] != j]
                        working_df = working_df[working_df['subject'] != j]
                        clustered.add(j)
    return None


f1 = collect_all_possibe_markers(sys.argv[2])
f2 = calculate_cov(f1, sys.argv[1])
build_connected_components(f2, sys.argv[4]).to_csv(sys.argv[5]+"bowtie2_cc.txt", sep="\t")

cluster(f1, f2, sys.argv[3], sys.argv[5]+"/")     
