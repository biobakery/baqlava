## this script is meant to be used manually, eg in a jupyter notebook where user checks that the genomes look as expected each time:

import pandas as pd
import os
import re

# Genome set: ICTV
ICTV = pd.read_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/VMR_MSL38_v1.csv")

ICTV = ICTV[ICTV['Genome coverage']!='No entry in Genbank']

ICTV = ICTV.copy()
ICTV['Virus REFSEQ accession'] = ICTV['Virus REFSEQ accession'].fillna(ICTV['Virus GENBANK accession'])

#adjust four lines for bad parsing of segmented genomes
ICTV_formatted = ICTV.replace("L: MN567049; MN567050; S:MN567048", "L: MN567049; L2: MN567050; S: MN567048")
ICTV_formatted = ICTV_formatted.replace("LK928904 (2253.0260)","LK928904")
ICTV_formatted = ICTV_formatted.replace("NC_003197","NC_010463")

#NC_003197 was replaced with record Felsduovirus Fels2, NC_010463 since ICTV record was incorrect (salmonella genome)

def split_out_all_accessions(df):
    df1 = df.copy()
    segmented = df1[df1['Virus REFSEQ accession'].str.contains(':')]
    
    segmented2 = segmented.copy()
    segmented2['Virus REFSEQ accession'] = segmented2['Virus REFSEQ accession'].apply(lambda x: re.split(r'; |, |;',x))
    segmented3 = segmented2.explode('Virus REFSEQ accession')
    segmented3 = segmented3[segmented3['Virus REFSEQ accession']!='']
    segmented3[['segmented','Virus REFSEQ accession']] = segmented3['Virus REFSEQ accession'].str.split(': |:', 1, expand=True)

    remainders = df1[~df1['Virus REFSEQ accession'].str.contains(':')]
    
    incomplete = remainders[remainders['Virus REFSEQ accession'].str.contains(';|,')]
    incomplete2 = incomplete.copy()
    incomplete2['Virus REFSEQ accession'] = incomplete2['Virus REFSEQ accession'].apply(lambda x: re.split(r'; |, |;|,',x))
    incomplete3 = incomplete2.explode('Virus REFSEQ accession')
    incomplete3['segmented'] = 'incomplete_genome_fragments'
    
    remainders2 = remainders[~remainders['Virus REFSEQ accession'].str.contains(';|,')]
    return pd.concat([segmented3, incomplete3, remainders2])

ICTV_formatted2 = split_out_all_accessions(ICTV_formatted)

#gather the actual genomes:
ICTV_accessions = list(ICTV_formatted2[['Virus REFSEQ accession']].drop_duplicates()['Virus REFSEQ accession'])

with open("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_V0_sandbox/20230815/esearch_ICTV_genomes.sh", "w") as bash:
    bash.write("#!/bin/bash")
    bash.write("\n")
    for i in ICTV_accessions:
        bash.write("esearch -db nucleotide -query " + i + " | efetch -format fasta >> /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/baqlava/database_generation/input/databases/ICTV.fasta")
        bash.write("\n")

ICTV_formatted2.to_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/baqlava/database_generation/input/reference_files/ICTV_ref.txt", sep="\t")

# actually just ran this in terminal but could also do it this way:
#os.system("bash /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_V0_sandbox/20230815/esearch_ICTV_genomes.sh")

# now open back up the collected ICTV fasta file and match headers with the accessions (they don't always match nicely unfortunately)

def format_reference(f1, f2, formatted):
    bash_accessions = []
    with open(f1, 'r') as bash:
        for i in bash:
            i = i.strip()
            #if len(i) > 0:
            i = i.split(" ")
            if len(i) == 1:
                pass
            else:
                bash_accessions.append(i[4])
    
    fasta_genomes = []
    fasta_accessions = []
    with open(f2, 'r') as fasta:
        for i in fasta:
            i = i.strip()
            if len(i) == 0:
                pass
            elif i[0]=='>':
                fasta_genomes.append(i[1:])
                fasta_accessions.append(i[1:].split(".")[0])
    
    df1 = pd.DataFrame({'accession':fasta_accessions, 'header':fasta_genomes})
    #df1['merge_type'] = 'fasta'
    df2 = pd.DataFrame({'accession':bash_accessions}) 
    #df2['merge_type'] = 'bash'
    df3 = pd.merge(df2, df1, on='accession', how='outer', indicator=True)
    # the code below was used to identify genomes that were not pulled with the bash script. 
    # These were obtained manually at https://www.ncbi.nlm.nih.gov/sites/batchentrez
    print('GENOMES TO ADD')
    for i in df3.query("_merge=='left_only'")['accession']:
        print(i)
    # the code below was used to identify genomes that were only in the downloaded genomes 
    # these were human and were removed manually from the database
    print('GENOMES TO REMOVE')
    for i in df3.query("_merge=='right_only'")['accession']:
        print(i)
    df4 = df3.query("_merge=='both'").rename(columns={'accession':'Virus REFSEQ accession'})
    
    df5 = pd.merge(df4[['Virus REFSEQ accession','header']], formatted.copy(), on='Virus REFSEQ accession', how='inner')
    print('All genomes in ref file match to an accession:', len(df5) == len(df4))
    # if the above does not print TRUE, manual adjustment may be needed. This will vary with each iteration of the ICTV VMR (metadata reference) 
    return df5
            

# first time prints out the genomes to add into the ICTV database (not pulled with bash)
format_reference("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_V0_sandbox/20230815/esearch_ICTV_genomes.sh", "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/baqlava/database_generation/input/databases/ICTV.fasta", ICTV_formatted2)
# second time all genomes should be present:
I1 = format_reference("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava_V0_sandbox/20230815/esearch_ICTV_genomes.sh", "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/baqlava/database_generation/input/databases/ICTV.fa", ICTV_formatted2)

I1.to_csv("/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlavav02/baqlava/database_generation/input/reference_files/ICTV_ref2.txt", sep="\t")
