import pandas as pd
import sys

#sysargv1 = reference_files/ICTV_ref2.txt
#sysargv2 = cluster/final_clusters.txt
#sysargv3 = final_databases/BAQLaVa_nucleotidedb_reference_file.txt
#sysargv4 = idmap

def assign_segmentation(f1):
    df1 = pd.read_csv(f1, sep="\t", index_col='Unnamed: 0')
    df2 = df1.dropna(subset=['segmented'])
    df3 = df2.groupby(['Species','Virus name(s)'], as_index=False).first()[['Species','Virus name(s)']].reset_index()
    df3['sg'] = df3['index'] + 1
    df3['sg'] = df3['sg'].astype(str)
    df3['segment_group'] = 'segment_group_' + df3['sg']
    df4 = df3[['Species','Virus name(s)','segment_group']]
    df5 = pd.merge(df4, df1, on=['Species','Virus name(s)'], how='outer')
    return df5

binomial_errors = pd.DataFrame({'n_observed_less_than_or_equal_to':[2,8,15,20,30,40,50,60,100], 'n_errors_permitted':[0,1,2,3,4,6,7,8,10]})

def format_reference_assign_taxonomy(ictv, clusters, ref, biner):
    df1 = pd.read_csv(clusters, sep="\t", index_col="Unnamed: 0")
    df2 = pd.read_csv(ref, sep="\t", index_col="Unnamed: 0")
    df3 = pd.merge(df2, ictv.copy(), left_on='Orig_Header', right_on='header', how='inner').drop_duplicates()
    df4 = pd.merge(df3, df1, left_on='New_Header', right_on='cluster_member', how='right').drop_duplicates()

    df5 = df4.copy()
    # account for a problematic ICTV genome that does not resolve at Realm 
    # remove this for new iterations of ICTV! But account for similar VGBs as necesssary:
    df5 = df5[df5['New_Header'] != 'BAQ00025852']
    
    # first get the easy VGBs: 
    # We only need to care about VGBs with an ICTV genome, for which we can even try to assign taxonomy:
    tax1 = df5.copy().query("database=='ICTV'")[['VGB']].drop_duplicates()
    # Next, assign taxonomy when there is consensus on species: 
    df6 = df5.copy().query("database=='ICTV'").groupby(['VGB','Species'], as_index=False).first().groupby('VGB', as_index=False).count().query("Species==1")[['VGB']]
    df7 = df5.copy().query("database=='ICTV'").groupby('VGB', as_index=False).first()
    tax2 = pd.merge(df6, df7, on='VGB', how='inner')[['VGB', 'Species', 'Realm','Kingdom', 
                                                      'Phylum', 'Class', 'Order', 'Family', 'Genus']]
    
    # now loop through the more difficult VGBs one by one:
    manloop = list(set(tax1['VGB']) - set(tax2['VGB']))
    tax_lookup = ictv.copy()[['Realm', 'Kingdom','Phylum', 'Class','Order', 'Family','Genus', 'Species']].drop_duplicates()
    
    VGBs = []
    
    Real_df = []
    King_df = []
    Phyl_df = []
    Cla_df = []
    Ord_df = []
    Fam_df = []
    Gen_df = []
    Spec_df = []
    
    
    for i in manloop:
        tempdf = df5.copy()[df5['VGB'] == i].query("database=='ICTV'")
        n_obs = len(tempdf)
        
        spec_lis = list(set(tempdf.dropna(subset=['Species'])['Species']))
        species = i + "_" + "_".join(spec_lis)
        species = species.replace(" ", "_")
        
        # first check to see if genus is concordant (most will be):
        # if not, use binomial error model
                
        if len(set(tempdf['Genus'])) == 1:
                    
            VGBs.append(i)
            Real_df.append(tempdf.iloc[0]['Realm']) 
            King_df.append(tempdf.iloc[0]['Kingdom']) 
            Phyl_df.append(tempdf.iloc[0]['Phylum']) 
            Cla_df.append(tempdf.iloc[0]['Class'])
            Ord_df.append(tempdf.iloc[0]['Order']) 
            Fam_df.append(tempdf.iloc[0]['Family']) 
            Gen_df.append(tempdf.iloc[0]['Genus'])
            Spec_df.append(species)
            
        else:
            ret = eval_genus_concordance(tempdf, biner, tax_lookup)
            if ret == "no concordance":
                ret = eval_family_concordance(tempdf, biner, tax_lookup)
                if ret == "no concordance":
                    ret = eval_order_concordance(tempdf, biner, tax_lookup)
                    if ret == "no concordance":
                        ret = eval_class_concordance(tempdf, biner, tax_lookup)
                        if ret == "no concordance":
                            print("ERROR: needs to resolve above class, code needs to be updated")
                            return "ERROR: needs to resolve above class, code needs to be updated"
                        else:
                            VGBs.append(i)
                            Real_df.append(ret[0]) 
                            King_df.append(ret[1]) 
                            Phyl_df.append(ret[2]) 
                            Cla_df.append(ret[3])
                            Ord_df.append("c__" + ret[3]+"_s__"+species) 
                            Fam_df.append("c__" + ret[3]+"_s__"+species) 
                            Gen_df.append("c__" + ret[3]+"_s__"+species)
                            Spec_df.append(species)
                            
                    else:
                        VGBs.append(i)
                        Real_df.append(ret[0]) 
                        King_df.append(ret[1]) 
                        Phyl_df.append(ret[2]) 
                        Cla_df.append(ret[3])
                        Ord_df.append(ret[4]) 
                        Fam_df.append("o__" + ret[4]+"_s__"+species) 
                        Gen_df.append("o__" + ret[4]+"_s__"+species)
                        Spec_df.append(species)
                        
                else:
                    VGBs.append(i)
                    Real_df.append(ret[0]) 
                    King_df.append(ret[1]) 
                    Phyl_df.append(ret[2]) 
                    Cla_df.append(ret[3])
                    Ord_df.append(ret[4]) 
                    Fam_df.append(ret[5]) 
                    Gen_df.append("f__" + ret[5]+"_s__"+species)
                    Spec_df.append(species)
            else:
                VGBs.append(i)
                Real_df.append(ret[0]) 
                King_df.append(ret[1]) 
                Phyl_df.append(ret[2]) 
                Cla_df.append(ret[3])
                Ord_df.append(ret[4]) 
                Fam_df.append(ret[5]) 
                Gen_df.append(ret[6])
                Spec_df.append(species)
                
            #print(ret)
            
    df8 = pd.DataFrame({'VGB':VGBs, 'Realm':Real_df, 'Kingdom':King_df, 
                        'Phylum':Phyl_df, 'Class':Cla_df, 'Order':Ord_df, 
                        'Family':Fam_df, 'Genus':Gen_df, 'Species':Spec_df})
    df9 = pd.concat([df8, tax2])
    df10 = pd.merge(df1[['VGB']].drop_duplicates(), df9, on='VGB', how='outer')
    df10['Species'] = df10['Species'].fillna(df10['VGB'])
    df10['Genus'] = df10['Genus'].fillna(df10['VGB'])
    return df5, df10
  
    
    
def eval_genus_concordance(df, biner, tax_lookup):
    df1 = df.copy()
    df2 = df1.dropna(subset=['Genus'])
    n_obs = len(df2)
    
    if n_obs == 0:
        return "no concordance"
    
    for j in range(len(biner)):
        if biner.iloc[j]['n_observed_less_than_or_equal_to'] >= n_obs:
            n_err_allowed = biner.iloc[j]['n_errors_permitted']
            #print(n_obs, 'genomes observed, ', n_err_allowed, 'errors allowed to call taxonomy')
            break
    
    # check to see if genus can meet genus assignment with binomial error model:
    n_genus_concord = df2.dropna(subset=['Genus']).groupby('Genus', as_index=False).count().max()['Species']
    n_genus_err = n_obs - n_genus_concord
    if int(n_genus_err) <= int(n_err_allowed):
                
        genus = df2[['Genus','Species']].groupby('Genus', as_index=False).count().sort_values(by='Species', ascending = False).iloc[0]['Genus']
        tax = tax_lookup.copy()[tax_lookup['Genus'] == genus]
        return [tax.iloc[0]['Realm'], tax.iloc[0]['Kingdom'], tax.iloc[0]['Phylum'], tax.iloc[0]['Class'], tax.iloc[0]['Order'], tax.iloc[0]['Family'], genus]
    
    else:
        
        return "no concordance"
    
def eval_family_concordance(df, biner, tax_lookup):
    df1 = df.copy()
    df2 = df1.dropna(subset=['Family'])
    n_obs = len(df2)
    
    if n_obs == 0:
        return "no concordance"
    
    for j in range(len(biner)):
        if biner.iloc[j]['n_observed_less_than_or_equal_to'] >= n_obs:
            n_err_allowed = biner.iloc[j]['n_errors_permitted']
            #print(n_obs, 'genomes observed, ', n_err_allowed, 'errors allowed to call taxonomy')
            break
    
    # check to see if genus can meet genus assignment with binomial error model:
    n_genus_concord = df2.groupby('Family', as_index=False).count().max()['Species']
    n_genus_err = n_obs - n_genus_concord
    if int(n_genus_err) <= int(n_err_allowed):
                
        family = df2[['Family','Species']].groupby('Family', as_index=False).count().sort_values(by='Species', ascending = False).iloc[0]['Family']
        tax = tax_lookup.copy()[tax_lookup['Family'] == family]
        return [tax.iloc[0]['Realm'], tax.iloc[0]['Kingdom'], tax.iloc[0]['Phylum'], tax.iloc[0]['Class'], tax.iloc[0]['Order'], family]
    
    else:
        
        return "no concordance"
    
def eval_order_concordance(df, biner, tax_lookup):
    df1 = df.copy()
    df2 = df1.dropna(subset=['Order'])
    n_obs = len(df2)
    
    if n_obs == 0:
        return "no concordance"
    
    for j in range(len(biner)):
        if biner.iloc[j]['n_observed_less_than_or_equal_to'] >= n_obs:
            n_err_allowed = biner.iloc[j]['n_errors_permitted']
            #print(n_obs, 'genomes observed, ', n_err_allowed, 'errors allowed to call taxonomy')
            break
    
    # check to see if genus can meet genus assignment with binomial error model:
    n_genus_concord = df2.groupby('Order', as_index=False).count().max()['Species']
    n_genus_err = n_obs - n_genus_concord
    if int(n_genus_err) <= int(n_err_allowed):
                
        order = df2[['Order','Species']].groupby('Order', as_index=False).count().sort_values(by='Order', ascending = False).iloc[0]['Order']
        tax = tax_lookup.copy()[tax_lookup['Order'] == order]
        return [tax.iloc[0]['Realm'], tax.iloc[0]['Kingdom'], tax.iloc[0]['Phylum'], tax.iloc[0]['Class'], order]
    
    else:
        
        return "no concordance"
    
def eval_class_concordance(df, biner, tax_lookup):
    df1 = df.copy()
    df2 = df1.dropna(subset=['Class'])
    n_obs = len(df2)
    
    if n_obs == 0:
        return "no concordance"
    
    for j in range(len(biner)):
        if biner.iloc[j]['n_observed_less_than_or_equal_to'] >= n_obs:
            n_err_allowed = biner.iloc[j]['n_errors_permitted']
            #print(n_obs, 'genomes observed, ', n_err_allowed, 'errors allowed to call taxonomy')
            break
    
    # check to see if genus can meet genus assignment with binomial error model:
    n_genus_concord = df2.groupby('Class', as_index=False).count().max()['Species']
    n_genus_err = n_obs - n_genus_concord
    if int(n_genus_err) <= int(n_err_allowed):
                
        Vclass = df2[['Class','Species']].groupby('Class', as_index=False).count().sort_values(by='Class', ascending = False).iloc[0]['Class']
        tax = tax_lookup.copy()[tax_lookup['Class'] == Vclass]
        return [tax.iloc[0]['Realm'], tax.iloc[0]['Kingdom'], tax.iloc[0]['Phylum'], Vclass]
    
    else:
        
        return "no concordance"


def account_for_segmentation_step1(r1, r2):
    df1 = r1[['VGB','Species','Virus name(s)','segment_group','segmented']].drop_duplicates().dropna(subset=['segment_group'])
    
    # ignore incomplete_genome_fragments for now: come back to later
    #df2 = df1.query("segmented!='incomplete_genome_fragments'")
    df3 = pd.merge(df1.rename(columns={'Species':'ictv_species'}), r2.copy(), on='VGB', how='left')[['VGB', 'ictv_species', 'Virus name(s)', 'segment_group','Realm', 'Kingdom','Phylum', 'Class','Order', 'Family', 'Genus', 'Species']].drop_duplicates()
    
    # len = 4534 : VGB - species group lines
    df4 = df3.copy()[['VGB','segment_group']].drop_duplicates().groupby('VGB', as_index=False).count().rename(columns={'segment_group':'n_segment_groups_VGB_in'})
    df5 = pd.merge(df3, df4, on='VGB', how='outer')
    df6 = df5.query("n_segment_groups_VGB_in>1")[['segment_group']].drop_duplicates()
    df7 = pd.merge(df5, df6, on='segment_group', how='outer', indicator=True).drop_duplicates()
    
    # this should be a handful of genomes to check by hand - likely they are not actually segmented
    # because if they _were_ segmented but with unique VGBs, 
    # they would have been named so that their VGB name and their ICTV name would match
    # what happened in most of these cases is that one of their representative genomes was entered 
    # in a way so that it _looks_ like a segmented genome, but it is not actually
    
    #df8 = df7.query("_merge=='left_only'")
    #df9 = df8[df8['ictv_species'] != df8['Species']]
    
    #for i in range(len(df9)):
    #    print(df9.iloc[i]['ictv_species'])
    #    print(df9.iloc[i]['Species'])
    #    print("\n")
    
    # now fix the df as needed based on manual investigation above:
    df10 = df7.copy().query("ictv_species=='Bracoviriform congregatae'")
    df10['Species'] = 'VGB_5417_Bracoviriform_glomeratae_Bracoviriform_congregatae_Bracoviriform_flavipedis_Bracoviriform_rubeculae_VGB_25131_Bracoviriform_congregatae_Bracoviriform_rubeculae_VGB_27158_Bracoviriform_flavipedis_Bracoviriform_congregatae_Bracoviriform_rubeculae_VGB_29878_Bracoviriform_rubeculae_Bracoviriform_marginiventris_Bracoviriform_congregatae_Bracoviriform_glomeratae_Bracoviriform_flavipedis'
    df11 = pd.concat([df10, df7.copy().query("ictv_species!='Bracoviriform congregatae'")])
    
    df12 = df11.copy().query("ictv_species=='Jatropha mosaic virus'")
    df12['Species'] = 'VGB_77223_Jatropha_mosaic_virus_Tobacco_mottle_leaf_curl_virus'
    df13 = pd.concat([df12, df11.copy().query("ictv_species!='Jatropha mosaic virus'")])
    
    df14 = df13.copy().query("ictv_species=='Tomato yellow leaf curl Thailand virus'")
    df14['Species'] = 'VGB_19033_Tomato_yellow_leaf_curl_Thailand_virus_Squash_leaf_curl_Yunnan_virus'
    df15 = pd.concat([df14, df13.copy().query("ictv_species!='Tomato yellow leaf curl Thailand virus'")])
    
    df15 = df15.query("ictv_species!='Cotton leaf curl Gezira virus'")
    df15 = df15.query("ictv_species!='Honeysuckle yellow vein virus'")
    df15 = df15.query("ictv_species!='Tomato leaf curl Sri Lanka virus'")
    df15 = df15.query("ictv_species!='Tomato leaf curl Arusha virus'")
    df15 = df15.query("ictv_species!='Eupatorium yellow vein virus'")
    df15 = df15.query("ictv_species!='Sida yellow vein Vietnam virus'")
    df15 = df15.query("ictv_species!='Tomato leaf curl Ghana virus'")
    df15 = df15.query("ictv_species!='Papaya leaf curl Guandong virus'")
    # final length 4526 lines
    
    return df15
    
def account_for_segmentation_step2(step1_df):
    df1 = step1_df.copy().drop_duplicates()
    # easy groups first: 
    df2 = df1.copy().query("_merge=='left_only'")
    # len = 3861
    #df2['segment_handling_group'] = 'group1'
    # group 1 segment groups can be simply averaged over
    df2['VGB_specificity'] = 'group1'
    df2['VGB_specificity_type'] = 'subspecies'
    # group 1 segment groups can be simply averaged over
    
    # more difficult groups:
    df3 = df1.copy().query("_merge=='both'")
    lis = list(set(df3['segment_group']))
    
    retdf = pd.DataFrame({})
    processed = []
    
    for i in lis:
        if i in processed:
            #print('passing:', i)
            pass
        else:
            tempdf = df3.copy()[df3['segment_group'] == i]
            tempdf2 = pd.merge(df3.copy(), tempdf[['VGB']].drop_duplicates(), on='VGB', how='inner')
            tempdf3 = pd.merge(df3.copy(), tempdf2[['segment_group']].drop_duplicates(), on='segment_group', how='inner')
            tempdf4 = pd.merge(df3.copy(), tempdf3[['VGB']].drop_duplicates(), on='VGB', how='inner')
            tempdf5 = pd.merge(df3.copy(), tempdf4[['segment_group']].drop_duplicates(), on='segment_group', how='inner')
            tempdf6 = pd.merge(df3.copy(), tempdf5[['VGB']].drop_duplicates(), on='VGB', how='inner')
            tempdf7 = pd.merge(df3.copy(), tempdf6[['segment_group']].drop_duplicates(), on='segment_group', how='inner')
            tempdf8 = pd.merge(df3.copy(), tempdf7[['VGB']].drop_duplicates(), on='VGB', how='inner')
            tempdf9 = pd.merge(df3.copy(), tempdf8[['segment_group']].drop_duplicates(), on='segment_group', how='inner')
            tempdf10 = pd.merge(df3.copy(), tempdf9[['VGB']].drop_duplicates(), on='VGB', how='inner')
            tempdf11 = pd.merge(df3.copy(), tempdf10[['segment_group']].drop_duplicates(), on='segment_group', how='inner')
            tempdf12 = pd.merge(df3.copy(), tempdf11[['VGB']].drop_duplicates(), on='VGB', how='inner')
            if len(tempdf12) != len(tempdf10):
                print('ERROR')
                break
            tempdf12 = tempdf12.drop(columns={"_merge"}).fillna('unannotated')#[['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus']].fillna('unannotated')
            # this set of VGBs will now be a new segmentation group, 
            # and we are just updating the annotation to be helpful in assigning taxonomy
            # gather information:
            i1 = len(set(tempdf12['segment_group']))
            i2 = len(set(tempdf12['ictv_species']))
            if len(set(tempdf12['Genus'])) > 1:
                print('FUCK')
                
            if tempdf12.max()['n_segment_groups_VGB_in'] == i1:
                # we have at least one VGB that is conserved in all of the segment groups,
                # which will serve as a conserved marker,
                # and possibly some that are more unique to call taxonomy from
                if i2 == 1:
                    # all segment groups are characterized to the same ICTV species, so
                    # we can potentially differentiate at the strain / subspecies level
                    bdf1 = tempdf12.copy()[tempdf12['n_segment_groups_VGB_in'] == i1][['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus','ictv_species']].drop_duplicates()
                    bdf2 = tempdf12.copy()[tempdf12['n_segment_groups_VGB_in'] != i1].groupby(['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus','ictv_species'])['Virus name(s)'].apply('_'.join).reset_index()
                    bdf3 = pd.concat([bdf1, bdf2])
                    bdf3['Virus name(s)'] = bdf3['Virus name(s)'].fillna('conserved')
                    bdf3 = bdf3.rename(columns={'ictv_species':'Species','Virus name(s)':'VGB_specificity'})
                    bdf3['VGB_specificity_type'] = 'subspecies'
                    #bdf3['segment_handling_group'] = 'group2'
                    bdf3['segment_group'] = tempdf12.iloc[0]['segment_group']
                    processed += list(set(tempdf12['segment_group']))
                    retdf = pd.concat([retdf, bdf3])
                    #print(set(bdf3['Species']))
                    #print(set(tempdf12['Species']))
                    #print(bdf3.columns)
    
                else:
                    # segment groups are characterized to different ICTV species, so
                    # we can potentially differentiate at species level but not strain / subspecies level
                    bdf1 = tempdf12.copy()[tempdf12['n_segment_groups_VGB_in'] == i1].groupby(['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus'])['ictv_species'].apply('_'.join).reset_index()
                    bdf1['VGB_specificity'] = 'conserved'
                    bdf1 = bdf1.rename(columns={'ictv_species':'Species'})
                    bdf2 = tempdf12.copy()[tempdf12['n_segment_groups_VGB_in'] != i1].groupby(['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus'])['ictv_species'].apply('_'.join).reset_index()
                    bdf2 = bdf2.rename(columns={'ictv_species':'VGB_specificity'})
                    bdf3 = pd.concat([bdf1, bdf2])
                    bdf3['Species'] = bdf3['Species'].fillna(bdf1.iloc[0]['Species'])
                    bdf3['segment_group'] = tempdf12.iloc[0]['segment_group']
                    bdf3['VGB_specificity_type'] = 'species'
                    processed += list(set(tempdf12['segment_group']))
                    retdf = pd.concat([retdf, bdf3])
                    #print(set(bdf3['Species']))
                    #print(set(tempdf12['Species']))
                    #print(bdf3.columns)
            else:
                # we do not have any VGBs that are conserved in all of the segment groups
                if i2 > 1:
                    tdf1 = tempdf12[['VGB','ictv_species']].drop_duplicates()
                    tdf2 = tdf1.groupby('VGB', as_index=False).count().query("ictv_species==1")[['VGB']]
                    bdf1 = pd.merge(tdf2.copy(), tempdf12.copy(), on='VGB', how='inner')[['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus','Species']].drop_duplicates().rename(columns={'Species':'VGB_specificity'})
                    bdf2 = pd.merge(tdf2.copy(), tempdf12.copy(), on='VGB', how='outer', indicator=True).query("_merge=='right_only'")[['ictv_species','VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus']].drop_duplicates().groupby(['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus'])['ictv_species'].apply('_'.join).reset_index().rename(columns={'ictv_species':'VGB_specificity'})
                    bdf3 = pd.concat([bdf1, bdf2])
                    spec_list = '_'.join(list(set(tempdf12['ictv_species'])))
                    bdf3['Species'] = spec_list
                    bdf3['segment_group'] = tempdf12.iloc[0]['segment_group']
                    bdf3['VGB_specificity_type'] = 'species'
                    processed += list(set(tempdf12['segment_group']))
                    retdf = pd.concat([retdf, bdf3])
                else:
                    #print(i1, i2)
                    print('many species groups but one species')
    retdf = pd.concat([retdf, df2[['VGB','Realm', 'Kingdom','Phylum', 'Class','Order','Family','Genus','Species','segment_group', 'VGB_specificity', 'VGB_specificity_type', 'Virus name(s)']]])
    return retdf.reset_index()[['VGB', 'Realm', 'Kingdom','Phylum', 'Class','Order','Family', 'Genus', 'Species', 'VGB_specificity', 'segment_group', 'VGB_specificity_type','Virus name(s)']]
                
def combine_taxonomy_and_segmentation_and_genomes_for_markers(r1, r2, seg, idm, r3):
    df1 = r1.copy()[['VGB','cluster_member']].drop_duplicates()
    df2 = r2.copy()
    df3 = seg.copy()[['VGB']].drop_duplicates()
    df4 = pd.merge(df2, df3, on='VGB', how='outer', indicator=True).query("_merge=='left_only'")
    df5 = pd.concat([df4, seg.copy()])
    df6 = pd.merge(df1, df5, on='VGB', how='inner')
    
    df7 = pd.read_csv(idm, sep="\t", names=['marker','marker1','marker_len'])
    clus = []
    scs = []
    for i in df7['marker']:
        clus.append(i.split("_")[0])
        if "_sc" in i:
            scs.append('yes')
        else:
            scs.append('no')
    df7['cluster_member'] = clus
    df7['supercluster_marker'] = scs
    
    df8 = pd.merge(df6.drop(columns={"_merge"}), df7[['cluster_member','marker','marker_len','supercluster_marker']], on='cluster_member', how='outer')
    df8['marker'] = df8['marker'].fillna('no_marker')
    
    df9 = pd.read_csv(r3, sep="\t", index_col="Unnamed: 0").groupby(["VGB","cluster_rep"], as_index=False).first().groupby(["VGB"], as_index=False).count().query("cluster_rep>1")
    df9['supercluster_VGB'] = 'yes'
    
    df10 = pd.merge(df8, df9[["VGB", "supercluster_VGB"]], on="VGB", how="outer")
    df10['supercluster_VGB'] = df10['supercluster_VGB'].fillna('no')
    
    return df10                 

ictv1 = assign_segmentation(sys.argv[1])
ref1, ref2 = format_reference_assign_taxonomy(ictv1, sys.argv[2], sys.argv[3], binomial_errors)
s1 = account_for_segmentation_step1(ref1, ref2)
s2 = account_for_segmentation_step2(s1)
combine_taxonomy_and_segmentation_and_genomes_for_markers(ref1, ref2, s2, sys.argv[4], sys.argv[2]).to_csv(sys.argv[5], sep="\t")
