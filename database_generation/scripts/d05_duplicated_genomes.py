import pandas as pd
import sys

#argv1: sam alignments to self
#argv2: location to save file 

def process_duplicates(file):
    d1 = pd.read_csv(file, sep="\t")
    d2 = d1.groupby("query", as_index=False).count()
    d2[['query','read']] = d2['query'].str.split("_", expand=True)
    d3 = d2.groupby("query", as_index=False).mean()
    d4 = d3[d3['subject']>=1.5]
    return d4[['query']].rename(columns={'query':'genome'})
    
process_duplicates(sys.argv[1]).to_csv(sys.argv[2] +"duplicated_genomes.txt", sep="\t")
