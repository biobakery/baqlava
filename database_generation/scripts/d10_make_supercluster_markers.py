import pandas as pd
import sys

#sysargv1: markers/supercluster_100mers.txt"
#sysargv2: cluster/final_clusters.txt"
#sysargv3: bowtie2/BAQLaVa_dereplicated_100mers.fasta"
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
        
        marker_lis = markers_dict[i]
        marker_lis.sort()
        
        for j in marker_lis:
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
        with open(loc+"supercluster_markers.fasta", 'w') as out:
            writer = 'off'
            cur_marker = ''
            for i in all_reads:
                i = i.strip()
                if len(i) == 0:
                    pass
                elif marker_counter == len(df2):
                    # basically we have gone through every 100mer, so now if we try to index
                    # we will not be in an .iloc[] that is within the dataframe
                    return None
                elif i[0] == '>':
                    if i[1:] == df2.iloc[marker_counter]['hundred_mer']:
                        if cur_marker != df2.iloc[marker_counter]['marker_number']: 
                            #print("match hit", df2.iloc[marker_counter])
                            cur_marker = df2.iloc[marker_counter]['marker_number']
                            writer = 'on'
                            out.write("\n")
                            out.write("\n")
                            out.write(i.split("_")[0] + "_" + str(df2.iloc[marker_counter]['markers']))
                            out.write("_sc")
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
    return None

def get_supercluster_genomes(file):
    df1 = pd.read_csv(file, sep="\t", index_col='Unnamed: 0')
    cc_size = df1.copy().groupby(['connected_component','cluster_rep'], as_index=False).first().groupby(['connected_component'], as_index=False).count().query("cluster_rep>1")
    df2 = df1[['cluster_rep','connected_component']].drop_duplicates()
    df3 = pd.merge(df2, cc_size, on='connected_component', how='inner')
    return df3.rename(columns={'cluster_rep_x':'genome'})[['genome','connected_component']]

def collect_marker_stats(mdf_2, gen_df):
    df1 = mdf_2.copy()
    genome = []
    marker_len = []
    n_markers = []
    df2 = df1.groupby(['genome','markers'], as_index=False).count().rename(columns={'hundred_mer':'marker_len'})
    df3 = pd.merge(df2, gen_df, on='genome', how='inner')
    
    df3_nomin = df3.copy().groupby('connected_component').sum().rename(columns={'marker_len':'2hbp_addit_marker_len'})
    df3_300 = df3.copy().query("marker_len>=3").groupby('connected_component').sum().rename(columns={'marker_len':'3hbp_addit_marker_len'})
    df3_400 = df3.copy().query("marker_len>=4").groupby('connected_component').sum().rename(columns={'marker_len':'4hbp_addit_marker_len'})
    df3_500 = df3.copy().query("marker_len>=5").groupby('connected_component').sum().rename(columns={'marker_len':'5hbp_addit_marker_len'})
    df3_600 = df3.copy().query("marker_len>=6").groupby('connected_component').sum().rename(columns={'marker_len':'6hbp_addit_marker_len'})
    df3_700 = df3.copy().query("marker_len>=7").groupby('connected_component').sum().rename(columns={'marker_len':'7hbp_addit_marker_len'})
    df3_800 = df3.copy().query("marker_len>=8").groupby('connected_component').sum().rename(columns={'marker_len':'8hbp_addit_marker_len'})
    df3_900 = df3.copy().query("marker_len>=9").groupby('connected_component').sum().rename(columns={'marker_len':'9hbp_addit_marker_len'})
    df3_1000 = df3.copy().query("marker_len>=10").groupby('connected_component').sum().rename(columns={'marker_len':'1kbp_addit_marker_len'})
    df3_1500 = df3.copy().query("marker_len>=15").groupby('connected_component').sum().rename(columns={'marker_len':'1p5kbp_addit_marker_len'})
    df3_2000 = df3.copy().query("marker_len>=20").groupby('connected_component').sum().rename(columns={'marker_len':'2kbp_addit_marker_len'})
    
    return pd.concat([df3.groupby('connected_component').first(), df3_nomin, df3_300, df3_400, df3_500, df3_600, df3_700, df3_800, df3_900, df3_1000, df3_1500, df3_2000], axis=1).fillna(0)


marker_dct = make_dictionary(sys.argv[1])
markerdf1 = get_contiguous_markers_faster(marker_dct)
markerdf2 = format_markers_faster(markerdf1)
markerdf2.groupby(['genome','markers'], as_index=False).count().rename(columns={'hundred_mer':'length_in_100bp'}).to_csv(sys.argv[4] + "/supercluster_marker_lengths.txt", sep="\t")
write_markers_faster(markerdf2, sys.argv[3], sys.argv[4]+"/")
genomes1 = get_supercluster_genomes(sys.argv[2])
collect_marker_stats(markerdf2, genomes1).to_csv(sys.argv[4]+"/supercluster_marker_stats.txt")

