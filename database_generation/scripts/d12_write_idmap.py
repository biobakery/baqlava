import pandas as pd
import sys

#sysargv1: markers/clusterrep_marker_lengths.txt
#sysargv2: markers/supercluster_marker_lengths.txt
#sysargv4: markers/all_markers.fasta

def write_idmap(f1, f2, f3):
    df1 = pd.read_csv(f1, sep="\t", index_col="Unnamed: 0")
    df1['markers'] = df1['markers'].astype(str)
    df1['c1'] = df1['genome'] + "_" + df1['markers']
    
    df2 = pd.read_csv(f2, sep="\t", index_col="Unnamed: 0")
    df2['markers'] = df2['markers'].astype(str)
    df2['c1'] = df2['genome'] + "_" + df2['markers'] + "_sc"
    
    df3 = pd.concat([df1,df2])
    df3['c2'] = df3['c1']
    df3['c3'] = df3['length_in_100bp'] * 100
    
    df4 = df3[['c1','c2','c3']]
    
    final_marker_set = []
    with open(f3, "r") as fa:
        for i in fa:
            i = i.strip()
            if len(i) == 0:
                pass
            elif i[0] == ">":
                final_marker_set.append(i[1:])
            else:
                pass
    
    df5 = pd.merge(pd.DataFrame({'c1':final_marker_set}), df4, on='c1', how='inner')
    
    return df5


write_idmap(sys.argv[1], sys.argv[2], sys.argv[3]).to_csv(sys.argv[4], sep="\t", header=False, index=False)
