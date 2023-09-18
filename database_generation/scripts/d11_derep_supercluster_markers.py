import pandas as pd
import os
import sys

#sysargv1 = /markers/supercluster_markers.fasta
#sysargv2 = /markers/SC_clu.tsv

def generate_database(arg1, arg2, arg3):
    df1 = pd.read_csv(arg2, sep="\t", names = ["cluster_rep", "cluster_member"])
    lis1 = list(df1[["cluster_rep"]].drop_duplicates()["cluster_rep"])
    keeper = 0
    genome = ""
    string = ""
    with open(arg1, "r") as all_genomes:
        with open(arg3, "w") as write_db:
            for i in all_genomes:
                i = i.strip()
                if len(i) == 0:
                    pass
                elif i[0] == '>':
                    # write out previous genome 
                    if keeper == 0:
                        pass
                    elif keeper == 1:
                        write_db.write(">" + genome)
                        write_db.write("\n")
                        write_db.write(string)
                        write_db.write("\n")
                    # reset values for current genome & header 
                    keeper = 0
                    string = ""
                    genome = ""
                    # process header
                    i = i.split(">")[1]
                    if i in lis1:
                        genome = i
                        keeper = 1
                    else:
                        genome = ""
                        keeper = 0
                else:
                    if keeper == 0:
                        pass
                    elif keeper == 1:
                        string += i
                    
    return None

generate_database(sys.argv[1], sys.argv[2], sys.argv[3])
