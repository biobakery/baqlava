import pandas as pd
import os
import sys
import re

def compile_dict(file, counter_start, path):
    string = "BAQ"
    counter = counter_start
    header_orig = []
    header_new = []
    ATGC_base_percent = []
    GC_content = []
    length = []
    N_percent = []
    N_longest = []
    genome = ""
    with open(path + "BAQLaVa_intermediate.fa", "a") as joined_fasta: 
        with open(file, "r") as fasta:
            for i in fasta:
                i = i.strip()
                if len(i) == 0:
                    pass
                elif i[0] == ">":
                    # first use info from last run though to calculate stats:
                    try:
                        N_percent.append(genome.count("N")/len(genome))
                        NL = re.findall( "(N+)", genome )
                        if len(NL) > 0:
                            N_longest.append(len(max(NL)))
                        else:
                            N_longest.append(0)
                        count = genome.count("A") + genome.count("T") + genome.count("C") + genome.count("G") 
                        ATGC_base_percent.append(count/len(genome))
                        GC_content.append((genome.count("C") + genome.count("G"))/len(genome))
                        length.append(len(genome))
                    except:
                        if counter == counter_start:
                            pass
                        else:
                            print("ERROR: skipping unintended line")
                    # next, reset stats for this round: 
                    counter += 1
                    genome = ""
                    # finally, process this new header:
                    header_orig.append(i[1:])
                    header_new.append(string + str(counter).zfill(8))
                    joined_fasta.write(">" + string + str(counter).zfill(8))
                    joined_fasta.write("\n")
                else:
                    genome += i
                    joined_fasta.write(i)
                    joined_fasta.write("\n")
            #calculate stats from final round
            N_percent.append(genome.count("N")/len(genome))
            NL = re.findall( "(N+)", genome )
            if len(NL) > 0:
                N_longest.append(len(max(NL)))
            else:
                N_longest.append(0)
            count = genome.count("A") + genome.count("T") + genome.count("C") + genome.count("G") 
            ATGC_base_percent.append(count/len(genome))
            GC_content.append((genome.count("C") + genome.count("G"))/len(genome))
            length.append(len(genome))
        retdf = pd.DataFrame({'Orig_Header':header_orig,'New_Header':header_new, 'ATCG_bases':ATGC_base_percent, 'GC_content':GC_content, 'length':length, 'N_percent':N_percent, 'N_longest':N_longest})
        retdf['database'] = os.path.split(file)[1].split(".")[0]
    return retdf, counter

def loop_through_files(loc, path):
    counter = 0
    retdf = pd.DataFrame({'Orig_Header':[],'New_Header':[], 'ATCG_bases':[], 'GC_content':[], 'length':[], 'Ns':[]})
    for i in os.listdir(loc):
        x, y = compile_dict(loc + i, counter, path)
        print(i, " completed")
        print("current counter: ", counter)
        counter = y
        retdf = pd.concat([retdf, x])
    return retdf

def use_reffile_to_discard_genomes(df):
    df1 = df.copy()
    keep1 = df1.copy()[df1['database'] == 'ICTV']
    df2 = df1.copy()[df1['database'] != 'ICTV']
    df3 = df2[df2['N_percent']<0.05]
    df4 = df3[df3['N_longest']<10]
    df5 = df4[df4['ATCG_bases']>0.75]
    return pd.concat([keep1, df5])

# use df without problematic contigs to write final fasta file
def write_out_fasta_from_list(df, fasta, path):
    keeps = list(df['New_Header'])
    with open(fasta, "r") as fasta1:
        with open(sys.argv[3], "w") as fasta2:
            writer = 'off'
            for i in fasta1:
                i = i.strip()
                if len(i) == 0:
                    pass
                elif i[0] == '>':
                    writer = 'off'
                    if i[1:] in keeps:
                        writer = 'on'
                        fasta2.write(i)
                        fasta2.write("\n")
                    else:
                        writer = 'off'
                else:
                    if writer == 'on':
                        fasta2.write(i)
                        fasta2.write("\n")
                    else:
                        pass
    return None    

print(sys.argv[1], sys.argv[2])
genomes_1 = loop_through_files(sys.argv[1]+"/", sys.argv[2]+"/")
genomes_df = use_reffile_to_discard_genomes(genomes_1)
genomes_df.to_csv(sys.argv[2] + "/BAQLaVa_nucleotidedb_reference_file.txt", sep="\t")

print("writing final formatted database")
write_out_fasta_from_list(genomes_df, sys.argv[2] + "/BAQLaVa_intermediate.fa", sys.argv[2]+"/")

os.system("rm " + sys.argv[2] + "/BAQLaVa_intermediate.fa")
print("formatted database complete")
