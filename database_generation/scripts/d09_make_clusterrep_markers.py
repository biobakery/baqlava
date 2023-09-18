import pandas as pd
import sys

#sysargv1: markers/cluster_rep_candidate100mers.txt"
#sysargv2: bowtie2/BAQLaVa_dereplicated_100mers.fasta"
#sysargv3: final_clusters.txt"
#sysargv4: saving location

def make_dictionary(markers_file):
    df1 = pd.read_csv(markers_file, sep="\t", index_col='Unnamed: 0')
    
    mer_dct = {}
    for i in df1[['genome']].drop_duplicates()['genome']:
        mer_dct[i] = []
    for j in df1['read']:   
        mer_dct[j.split("_")[0]].append(int(j.split("_")[1]))
    
    return mer_dct


def get_contiguous_markers_faster(markers_dict):

    genomes = []
    marker_lists = []
    counter = 0
    
    for i in markers_dict:
        counter += 1
        #print("processing", i)
        markers = []
        marker_holder = []
        n = 0
        
        for j in markers_dict[i]:
            if n == 0:
                n = j
                marker_holder.append(j)
            else:
                if j == n + 1:
                    n = j
                    marker_holder.append(j)
                elif j >= n + 2:
                    if len(marker_holder) <= 1:
                        marker_holder = []
                        n = j
                        marker_holder.append(j)
                    else:
                        markers.append(marker_holder)
                        marker_holder = []
                        n = j
                        marker_holder.append(j)
        # end the loop:
        if len(marker_holder) >= 2:
            markers.append(marker_holder)
            marker_holder = []
            n = 0
        genomes.append(i)
        marker_lists.append(markers)
        #if len(markers) == 0:
        #    print('problematic genome:', i)
    return pd.DataFrame({'genome':genomes, 'markers':marker_lists})

def format_markers_faster(marker_df):
    
    df1 = marker_df.explode('markers').groupby('genome', as_index=False).count()
    df2 = marker_df.explode('markers').explode('markers').dropna()
    df2['markers'] = df2['markers'].astype(str)
    df2['hundred_mer'] = df2['genome'] + "_" + df2['markers']
    
    # for the 120k which just have 1 marker, move forward without any dataframe looping:
    df3 = df1.query("markers==1")
    df4 = pd.merge(df3, df2.copy()[['genome','hundred_mer']], on='genome', how='inner')
    
    # for the remaining 30k which have >1 marker, now we have to loop:
    df5 = df1.query("markers>1").sort_values(by='genome')
    df6 = pd.merge(df5[['genome']], marker_df.copy().explode('markers'), on='genome', how='inner')
    
    genomes = []
    markers = []
    ns = []
    
    for i in range(len(df5)):
        gen = df5.iloc[i]['genome']
        n = int(df5.iloc[i]['markers'])
        for j in range(n):
            genomes.append(gen)
            markers.append(gen + "_" + str(j+1))
            ns.append(j+1)
    
    df7 = pd.DataFrame({'genome_x':genomes, 'marker_no_1':markers, 'marker_no_2':ns})
    df8 = pd.concat([df6, df7], axis=1)
    
    # check that concat worked correctly:
    df9 = df8.copy()[df8['genome'] == df7['genome_x']]
    if len(df8) != len(df9):
        print('ERROR')
        return None
    
    #return df8
    df10 = df8.rename(columns={'markers':'hundred_mer_1'}).rename(columns={'marker_no_2':'markers'}).explode('hundred_mer_1')
    df10['hundred_mer_1'] = df10['hundred_mer_1'].astype(str)
    df10['hundred_mer'] = df10['genome'] + "_" + df10['hundred_mer_1']
    df11 = pd.concat([df10[['genome','markers','hundred_mer']], df4])
    df11['markers'] = df11['markers'].astype(int)
    
    print(len(df2.dropna())==len(df11))
    
    return df11.sort_values(by=['genome','markers'])
    
def write_markers_faster(df, file, loc):
    df1 = df.copy()
    df1['markers'] = df1['markers'].astype(str)
    df1['marker_number'] = df1['genome'] + "_" + df1['markers']

    df1['markers'] = df1['markers'].astype(int)
    df2 = df1.sort_values(by=['genome','markers']).reset_index(drop=True)
    
    marker_counter = 0
    
    with open(file, "r") as all_reads:
        with open(loc + "clusterrep_markers.fasta", 'w') as out:
            writer = 'off'
            cur_marker = ''
            for i in all_reads:
                i = i.strip()
                if len(i) == 0:
                    pass
                elif i[0] == '>':
                    if i[1:] == df2.iloc[marker_counter]['hundred_mer']:
                        if cur_marker != df2.iloc[marker_counter]['marker_number']: 
                            print("match hit", df2.iloc[marker_counter])
                            cur_marker = df2.iloc[marker_counter]['marker_number']
                            writer = 'on'
                            out.write("\n")
                            out.write("\n")
                            out.write(i.split("_")[0] + "_" + str(df2.iloc[marker_counter]['markers']))
                            out.write("\n")
                        marker_counter +=1
                    else:
                        writer = 'off'
                        cur_marker = ''
                else:
                    if writer == 'off':
                        pass
                    elif writer == 'on':
                        out.write(i)
            out.write("\n")
    return None



def get_cluster_genomes(file):
    df1 = pd.read_csv(file, sep="\t", index_col='Unnamed: 0')
    return df1[['cluster_rep']].drop_duplicates().rename(columns={'cluster_rep':'genome'})

def collect_marker_stats(mdf_2, gen_df):
    df1 = mdf_2.copy()
    genome = []
    marker_len = []
    n_markers = []
    df2 = df1.groupby(['genome','markers'], as_index=False).count().rename(columns={'hundred_mer':'marker_len'})
    df2_nomin = df2.copy().groupby('genome').sum().rename(columns={'marker_len':'2hbp_addit_marker_len'})
    df2_300 = df2.copy().query("marker_len>=3").groupby('genome').sum().rename(columns={'marker_len':'3hbp_addit_marker_len'})
    df2_400 = df2.copy().query("marker_len>=4").groupby('genome').sum().rename(columns={'marker_len':'4hbp_addit_marker_len'})
    df2_500 = df2.copy().query("marker_len>=5").groupby('genome').sum().rename(columns={'marker_len':'5hbp_addit_marker_len'})
    df2_600 = df2.copy().query("marker_len>=6").groupby('genome').sum().rename(columns={'marker_len':'6hbp_addit_marker_len'})
    df2_700 = df2.copy().query("marker_len>=7").groupby('genome').sum().rename(columns={'marker_len':'7hbp_addit_marker_len'})
    df2_800 = df2.copy().query("marker_len>=8").groupby('genome').sum().rename(columns={'marker_len':'8hbp_addit_marker_len'})
    df2_900 = df2.copy().query("marker_len>=9").groupby('genome').sum().rename(columns={'marker_len':'9hbp_addit_marker_len'})
    df2_1000 = df2.copy().query("marker_len>=10").groupby('genome').sum().rename(columns={'marker_len':'1kbp_addit_marker_len'})
    df2_1500 = df2.copy().query("marker_len>=15").groupby('genome').sum().rename(columns={'marker_len':'1p5kbp_addit_marker_len'})
    df2_2000 = df2.copy().query("marker_len>=20").groupby('genome').sum().rename(columns={'marker_len':'2kbp_addit_marker_len'})
    
    return pd.concat([gen_df.copy().groupby('genome').first(), df2_nomin, df2_300, df2_400, df2_500, df2_600, df2_700, df2_800, df2_900, df2_1000, df2_1500, df2_2000], axis=1).fillna(0)




marker_dct = make_dictionary(sys.argv[1])
markerdf1 = get_contiguous_markers_faster(marker_dct)
markerdf2 = format_markers_faster(markerdf1)
markerdf2.groupby(['genome','markers'], as_index=False).count().rename(columns={'hundred_mer':'length_in_100bp'}).to_csv(sys.argv[4] + "/clusterrep_marker_lengths.txt", sep="\t")
write_markers_faster(markerdf2, sys.argv[2], sys.argv[4]+"/")
genomes = get_cluster_genomes(sys.argv[3])
collect_marker_stats(markerdf2, genomes).to_csv(sys.argv[4]+"clusterrep_marker_stats.txt")
