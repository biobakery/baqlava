import pandas as pd
import sys

# sysargv1 : marker_information/marker_taxonomy.txt
# sysargv2 : markers/all_markers.fasta
# sysargv3 : markers/all_markers_translated.fasta
# sysargv4 : (out) name of demo nucleotide fasta file
# sysargv5 : (out) name of demo translated fasta file

def make_marker_df(file):
    df1 = pd.read_csv(file, sep="\t", index_col='Unnamed: 0')
    df2 = pd.concat([df1[df1["Species"].str.contains("Emesvirus")][['marker']].drop_duplicates(), 
                     df1[df1["Species"].str.contains("Sinsheimervirus")][['marker']].drop_duplicates(), 
                     df1[df1["Species"].str.contains("Tunavirus")][['marker']].drop_duplicates(), 
                     df1[df1["Species"].str.contains("Punavirus")][['marker']].drop_duplicates(), 
                     df1[df1["Species"].str.contains("Tequatrovirus")][['marker']].drop_duplicates()])
    return df2


def make_demo_nuc(mark_df, fasta, fout):
    df1 = mark_df.copy()
    mark_dct = {}
    for i in df1['marker']:
        mark_df[i] = 1
    writer = 'off'
    with open(fasta, "r") as fa:
        with open(fout, "w") as out:
            for j in fa:
                j = j.strip()
                if len(j) == 0:
                    pass
                elif j[0] == '>':
                    if j[1:] in mark_df:
                        writer = 'on'
                        out.write(j)
                        out.write("\n")
                    else:
                        writer = 'off'
                elif writer == 'on':
                    out.write(j)
                    out.write("\n")
                elif writer == 'off':
                    pass
                
    return None

def make_demo_trans(mark_df, fasta, fout):
    df1 = mark_df.copy()
    mark_dct = {}
    for i in df1['marker']:
        mark_df[i] = 1
    writer = 'off'
    with open(fasta, "r") as fa:
        with open(fout, "w") as out:
            for j in fa:
                j = j.strip()
                if len(j) == 0:
                    pass
                elif j[0] == '>':
                    if j[1:-4] in mark_df:
                        writer = 'on'
                        out.write(j)
                        out.write("\n")
                    else:
                        writer = 'off'
                elif writer == 'on':
                    out.write(j)
                    out.write("\n")
                elif writer == 'off':
                    pass
                
    return None


markers = make_marker_df(sys.argv[1])
make_demo_nuc(markers, sys.argv[2], sys.argv[4])
make_demo_trans(markers, sys.argv[3], sys.argv[5])
