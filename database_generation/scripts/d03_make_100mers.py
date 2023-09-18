import pandas as pd
import os
import sys

#sys.argv1 = full databse
#sys.argv2 = cluster output from anadama

def generate_kmers(arg1, arg2, arg3):
    df1 = pd.read_csv(arg2, sep="\t", names = ["cluster_rep", "cluster_member"])
    lis1 = list(df1[["cluster_rep"]].drop_duplicates()["cluster_rep"])
    keeper = 0
    genome = ""
    string = ""
    with open(arg1, "r") as all_genomes:
        with open(arg3 + "BAQLaVa_dereplicated_100mers.fasta", "w") as write_kmers:
            for i in all_genomes:
                i = i.strip()
                if len(i) == 0:
                    pass
                elif i[0] == '>':
                    # generate kmers from previous genome 
                    if keeper == 0:
                        pass
                    elif keeper == 1:
                        offset = (len(string)%100)//2
                        for j in range(len(string)//100):
                            write_kmers.write(">" + genome + "_" + str(j+1))
                            write_kmers.write("\n")
                            write_kmers.write(string[(j*100+offset):(j*100+100+offset)])
                            write_kmers.write("\n")
                    # reset values for current genome & header 
                    keeper = 0
                    string = ""
                    genome = ""
                    # process header
                    i = i.split(">")[1]
                    if i in lis1:
                        print(i)
                        genome = i
                        keeper = 1
                    else:
                        #print(i)
                        #print("passing")
                        genome = ""
                        keeper = 0
                else:
                    if keeper == 0:
                        pass
                    elif keeper == 1:
                        string += i
                    
    return None

generate_kmers(sys.argv[1], sys.argv[2], sys.argv[3]+"/")
