import pandas as pd
import sys
import re

#sysargv1 = cluster/final_clusters.txt
#sysargv2 = bowtie2/BAQLaVa_dereplicated_100mers.fasta
#sysargv3 = bowtie2/BAQ_100mers_mapped_sorted.sam
#sysargv4 = save location


def get_cluster_genomes_dct(file):
    df1 = pd.read_csv(file, sep="\t", index_col='Unnamed: 0')
    genome_dct = {}
    cluster_reps = {}
    all_genomes = {}
    for j in range(len(df1)):
        genome_dct[ df1.iloc[j]['cluster_member'] ] = [ df1.iloc[j]['cluster_rep'], df1.iloc[j]['VGB'] ]
        cluster_reps[df1.iloc[j]['cluster_rep']] = 1
        all_genomes[df1.iloc[j]['cluster_member']] = 1
    return df1, genome_dct, cluster_reps, all_genomes

def collect_all_possibe_markers(file, cluster_df):
    keep_genomes = {}
    df1 = cluster_df.copy()[['cluster_rep']].drop_duplicates()
    for i in df1['cluster_rep']:
        keep_genomes[i] = 1
    counter = 0
    gens = []
    reads = []
    with open(file, "r") as fasta:
        for i in fasta:
            i = i.strip()
            if i[0] == ">":
                i2 = i.split("_")[0]
                if keep_genomes.get(i2[1:]) == 1:
                    gens.append(i2[1:])
                    reads.append(i[1:])
                else:
                    pass
    return pd.DataFrame({'genome':gens, 'read':reads})

def calculate_percent_identity(cigar_string, md_field):
    
    sam_cigar_match_mismatch_indel_identifiers=["M","=","X","I","D"]
    sam_cigar_add_to_reference_identifiers=["M","D","N","=","X"]
    sam_md_field_identifier="MD:Z:"
    
    """
    Calculate the percent identity using the cigar string and md field from the sam file
    Returns the percent identity and the alignment length
    """
    
    match_numbers=re.compile("\d+")
    match_non_numbers=re.compile("\D+")
    
    # find the sets of numbers and identifers from the cigar string
    cigar_numbers=match_numbers.findall(cigar_string)
    cigar_identifiers=match_non_numbers.findall(cigar_string)

    # find the index for all of the match/mismatch/insert/delete
    match_mismatch_indel_index = []
    reference_length_index = []
    for index, cigar_identifier in enumerate(cigar_identifiers):
        #print(cigar_identifier)
        if cigar_identifier in sam_cigar_match_mismatch_indel_identifiers: #config.sam_cigar_match_mismatch_indel_identifiers:
            match_mismatch_indel_index.append(index)
        if cigar_identifier in sam_cigar_add_to_reference_identifiers: #config.sam_cigar_add_to_reference_identifiers:
            reference_length_index.append(index)
    

    # get reference length
    try:
        reference_length = int(sum([float(cigar_numbers[index]) for index in reference_length_index]))
    except (IndexError, ValueError):
        reference_length = 0

    
    # identify the total number of match/mismatch/indel
    try:
        match_mismatch_indel_count=sum([float(cigar_numbers[index]) for index in match_mismatch_indel_index])
    except (IndexError, ValueError):
        match_mismatch_indel_count=0.0 

    # remove the tag from the md field
    md_field=md_field.split(sam_md_field_identifier)[-1] #(config.sam_md_field_identifier)[-1]
    
    
    # find the sets of numbers from the md field
    md_field_numbers=match_numbers.findall(md_field)
    
    # sum the md field numbers to get the total number of matches
    try:
        matches=sum(int(n) for n in md_field_numbers)
    except ValueError:
        matches=0.0
    
    percent_identity=0.0
    if match_mismatch_indel_count > 0.0:
        percent_identity = 100.0 * ( matches / ( match_mismatch_indel_count * 1.0 ) )
        
    return percent_identity, match_mismatch_indel_count, reference_length 

def loop_thu_sam(file, df_100mers, gen_dct, clusterrep_idx, genome_idx, loc):
    cluster_100mers = set(df_100mers['read'])
    supercluster_100mers = set()
    discard_100mers = set()
    counter = 0
    with open(file, "r") as sam:
        with open(loc + "supercluster_alignments.txt", "w") as superclus_100mers:
                
            superclus_100mers.write("query_genome")
            superclus_100mers.write("\t")
            superclus_100mers.write("query_read")
            superclus_100mers.write("\t")
            superclus_100mers.write("subject")
            superclus_100mers.write("\t")
            superclus_100mers.write("percent_identity")
            superclus_100mers.write("\n")
                
            for i in sam:
                i = i.strip()
                if i[0] == '@':
                    pass
                else:
                    i = i.split()
                    query = i[0].split("_")[0]
                    subject = i[2]
                    if query == subject:
                        # mapped to self, expected behavior, ignore
                        pass
                    else:
                        #counter += 1
                        if clusterrep_idx.get(query) == 1:
                            # 100mer comes from a genome we want to markerize (a cluster representative)
                            if genome_idx.get(subject) == 1:
                                
                                #if counter == 10000:
                                #    cluster2 = cluster_100mers - discard_100mers
                                #    supercluster2 =  supercluster_100mers - discard_100mers - cluster_100mers
                                #    return cluster2, supercluster2
 
                                # subject is a genome we care about (eg not discarded as duplicate, etc)
                                
                                query_clus = gen_dct.get(query)[0]
                                subject_clus = gen_dct.get(subject)[0]
                                query_cc = gen_dct.get(query)[1]
                                subject_cc = gen_dct.get(subject)[1]
                                
                                if query_clus == subject_clus:
                                    # query and subject are in the same cluster
                                    # keep the 100mer as a condender for cluster markers
                                    pass
                                elif query_cc == subject_cc:
                                    # query and subject are in the same connected component
                                    # keep the 100mer for supercluster markers if alignment is >=95% ID
                                    CS = i[5]
                                    for m in i[11:]:
                                        if re.match("MD:Z:", m):
                                            MD = m
                                            pid = calculate_percent_identity(CS, MD)[0]
                                            if pid >= 95:
                                                cluster_100mers.discard(i[0])
                                                supercluster_100mers.add(i[0])
                                                
                                                # write out supercluster-enriched 100mers:
                                                ### IMPORTANT: these reads are not used in the pipeline currently, 
                                                ### however they COULD be used to create supercluster markers made from 
                                                ### 100mers that map more often within a supercluser than other 100mers
                                                superclus_100mers.write(query)
                                                superclus_100mers.write("\t")
                                                superclus_100mers.write(i[0])
                                                superclus_100mers.write("\t")
                                                superclus_100mers.write(i[2])
                                                superclus_100mers.write("\t")
                                                superclus_100mers.write(str(pid))
                                                superclus_100mers.write("\n")
                                else:
                                    # query 100mer mapped outside of the connected component 
                                    # throw out the 100mer, its is nonunique
                                    cluster_100mers.discard(i[0])
                                    discard_100mers.add(i[0])
                            else:
                                # mapped to a genome we threw out for some reason, ignore the alignment
                                pass
                        else:
                            # query read comes from a genome we dont care about (dublicate, clustered, etc)
                            pass

    cluster2 = cluster_100mers - discard_100mers
    supercluster2 =  supercluster_100mers - discard_100mers - cluster_100mers
    return cluster2, supercluster2


def write_out(st):
    df1 = pd.DataFrame({'read':list(st)})
    df1[['genome','read_sort']] = df1['read'].str.split("_", expand=True)
    df1['read_sort'] = df1['read_sort'].astype(int)
    df2 = df1.sort_values(by=['genome','read_sort'])
    return df2



genomes, dct1, dct2, dct3  = get_cluster_genomes_dct(sys.argv[1])
all_100mers = collect_all_possibe_markers(sys.argv[2], genomes)
cluster_set, supercluster_set = loop_thu_sam(sys.argv[3], all_100mers, dct1, dct2, dct3, sys.argv[4]+ "/")

write_out(cluster_set).to_csv(sys.argv[4] + "cluster_candidate100mers.txt", sep="\t")
write_out(supercluster_set).to_csv(sys.argv[4] + "supercluster_candidate100mers.txt", sep="\t")
