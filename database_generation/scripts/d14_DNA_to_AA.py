import os
import sys
import pandas as pd

# sysargv1 : marker_dir + "all_markers.fasta"
# sysargv2 : marker_info_dir + "idmap1.txt"
# sysargv3 : (out) marker_dir + "all_markers_translated.fasta"
# sysargv4 : (out) marker_info_dir + "idmap2.txt"
# sysargv5 : (out) marker_info_dir + "translated_markers_conversion.txt"


gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
      'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W'}

basepairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def translate_frameshifted( sequence ):
    translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
    return translate

def reverse_complement( sequence ):
    reversed_sequence = (sequence[::-1])
    rc = ''.join([basepairs.get(reversed_sequence[i], 'N') for i in range(len(sequence))])
    return rc

def split_to_70s( string ):
    counter = 0
    newstring=[]
    tempstring=""
    for i in string:
        counter += 1
        tempstring += i
        if counter == 70:
            newstring.append(tempstring)
            tempstring=""
            counter = 0
    newstring.append(tempstring)
    return newstring

def fasta_to_AA( fin, fout ):
    #counter = 0
    with open(fin, "r") as data:
        with open(fout, "w") as trans:
            for i in data: 
                i = i.strip()
                if len(i) == 0:
                    pass
                elif i[0] == ">":
                    #counter += 1
                    i = i.strip()
                    name = i
                else:
                    seq = str(i)
                    #POSITIVE FRAME 1
                    trans.write(name + "_P1")
                    trans.write("\n")
                    write = split_to_70s( translate_frameshifted(seq[0:]) )  # first frame
                    for j in write:
                        trans.write(j)
                        trans.write("\n")
                    #POSITIVE FRAME 2
                    trans.write(name + "_P2")
                    trans.write("\n")
                    write = split_to_70s( translate_frameshifted(seq[1:]) )  # second frame
                    for j in write:
                        trans.write(j)
                        trans.write("\n")
                    #POSITIVE FRAME 3
                    trans.write(name + "_P3")
                    trans.write("\n")
                    write = split_to_70s( translate_frameshifted(seq[2:]) )  # third frame
                    for j in write:
                        trans.write(j)
                        trans.write("\n")
                    #NEGATIVE FRAME 1
                    trans.write(name + "_M1")
                    trans.write("\n")
                    write = split_to_70s( translate_frameshifted(reverse_complement(seq)) )  # negative first frame
                    #for j in write:
                    for j in write:
                        trans.write(j)
                        trans.write("\n")
                    #NEGATIVE FRAME 2
                    trans.write(name + "_M2")
                    trans.write("\n")
                    write = split_to_70s( translate_frameshifted(reverse_complement(seq[:len(seq)-1])) )  # negative second frame
                    for j in write:
                        trans.write(j)
                        trans.write("\n")
                    #NEGATIVE FRAME 3
                    trans.write(name + "_M3")
                    trans.write("\n")
                    write = split_to_70s( translate_frameshifted(reverse_complement(seq[:len(seq)-2])) )  # negative second frame
                    for j in write:
                        trans.write(j)
                        trans.write("\n")
                    #reset name and sequence
                    name = ""
                    seq = ""
    return None

def write_idmap(file):
    df1 = pd.read_csv(file, sep="\t", names=['id1','id2','len'])
    df2 = df1.copy()
    df2['id1'] = df2['id1'] + "_P1"
    df2['id2'] = df2['id2'] + "_P1"
    df3 = df1.copy()
    df3['id1'] = df3['id1'] + "_P2"
    df3['id2'] = df3['id2'] + "_P2"
    df4 = df1.copy()
    df4['id1'] = df4['id1'] + "_P3"
    df4['id2'] = df4['id2'] + "_P3"
    df5 = df1.copy()
    df5['id1'] = df5['id1'] + "_M1"
    df5['id2'] = df5['id2'] + "_M1"
    df6 = df1.copy()
    df6['id1'] = df6['id1'] + "_M2"
    df6['id2'] = df6['id2'] + "_M2"
    df7 = df1.copy()
    df7['id1'] = df7['id1'] + "_M3"
    df7['id2'] = df7['id2'] + "_M3"

    df8 = pd.concat([df1, df2, df3, df4, df5, df6, df7])
    return df8

def make_trans_marker_conversion(file):
    
    df1 = pd.read_csv(file, sep="\t", names=['translated_marker','nucleotide_marker','len'])
    df2 = df1.copy()
    print(len(df2))
    df2['translated_marker'] = df2['translated_marker'] + "_P1"
    df3 = df1.copy()
    df3['translated_marker'] = df3['translated_marker'] + "_P2"
    df4 = df1.copy()
    df4['translated_marker'] = df4['translated_marker'] + "_P3"
    df5 = df1.copy()
    df5['translated_marker'] = df5['translated_marker'] + "_M1"
    df6 = df1.copy()
    df6['translated_marker'] = df6['translated_marker'] + "_M2"
    df7 = df1.copy()
    df7['translated_marker'] = df7['translated_marker'] + "_M3"
    df8 = pd.concat([df2, df3, df4, df5, df6, df7])[['translated_marker','nucleotide_marker']]
    
    return df8






fasta_to_AA(sys.argv[1], sys.argv[3])
write_idmap(sys.argv[2]).to_csv(sys.argv[4], sep="\t", header=False, index=False)
make_trans_marker_conversion(sys.argv[2]).to_csv(sys.argv[5], sep="\t")
