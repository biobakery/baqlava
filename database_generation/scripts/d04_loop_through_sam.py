import pandas as pd
import sys
import re


def make_clu_dct(file):
    df1 = pd.read_csv(file, sep="\t", names=['cluster_rep','cluster_member'])
    retdct = {}
    for i in range(len(df1)):
        rep = df1.iloc[i]['cluster_rep']
        memb = df1.iloc[i]['cluster_member']
        retdct[memb] = rep
    return retdct

def loop_thu_sam(file, clusdct, loc):
    with open(file, "r") as sam:
        with open(loc + "mapped_reads.txt", "w") as mapd:
            with open(loc + "self_mapped.txt", "w") as self_mapd:
                with open(loc + "not_marker_candidates.txt", "w") as not_markers:
                    cur_dct = {}
                    cur_query = ""
                    reads_set = set()
                    self_mapd.write("query")
                    self_mapd.write("\t")
                    self_mapd.write("subject")
                    self_mapd.write("\n")
                    mapd.write("query")
                    mapd.write("\t")
                    mapd.write("subject")
                    mapd.write("\t")
                    mapd.write("n_markers_mapped")
                    mapd.write("\t")
                    mapd.write("markers_mapped")
                    mapd.write("\n")
                    counter = 0
                    for i in sam:
                        i = i.strip()
                        if i[0] == '@':
                            pass
                        else:
                            i = i.split()
                            query = i[0].split("_")[0]
                            read = i[0].split("_")[1]
                            subject = i[2]
                            if query == subject:
                                self_mapd.write(i[0])
                                self_mapd.write("\t")
                                self_mapd.write(i[2])
                                self_mapd.write("\n")
                            else:
                                if query != cur_query:
                                    # write out previous dct to file & reset dict 
                                    for j in cur_dct:
                                        mapd.write(cur_query)
                                        mapd.write("\t")
                                        mapd.write(j)
                                        mapd.write("\t")
                                        mapd.write(str(len(cur_dct[j])))
                                        mapd.write("\t")
                                        mapd.write(str(list(cur_dct[j])))
                                        mapd.write("\n")
                                    #for k in reads_set:
                                    #    not_markers.write(k)
                                    #    not_markers.write("\n")
                                    cur_dct = {}
                                    cur_query = query
                                    #reads_set = set()
                                # now continue with adding reads to dictionary in memory    
                                if clusdct.get(subject) == query:
                                    pass
                                else:
                                    if subject in cur_dct:
                                        cur_dct[subject].add(i[0])
                                    else:
                                        cur_dct[subject] = set()
                                        cur_dct[subject].add(i[0])
                                
                                #CS = i[5]
                                #for m in i[11:]:
                                #    if re.match("MD:Z:", m):
                                #        MD = m
                                #        pid = calculate_percent_identity(CS, MD)[0]
                                #        if pid >= 95:
                                #            reads_set.add(i[0])

                    #end loop by adding final dct to file
                    for j in cur_dct:
                        mapd.write(cur_query)
                        mapd.write("\t")
                        mapd.write(j)
                        mapd.write("\t")
                        mapd.write(str(len(cur_dct[j])))
                        mapd.write("\t")
                        mapd.write(str(list(cur_dct[j])))
                        mapd.write("\n")
    return None
                                    

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



cluster_dct = make_clu_dct(sys.argv[2])                   
loop_thu_sam(sys.argv[1], cluster_dct, sys.argv[3] + "/") 
