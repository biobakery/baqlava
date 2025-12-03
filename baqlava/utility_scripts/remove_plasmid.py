#!/usr/bin/env python3

import pandas as pd
import sys
import os

"""
Usage:
    python remove_plasmid.py <input_directory> <remove_VGBs> <output_directory>

Arguments:
    input_directory     Path to location of BAQLaVa profiles to remove plasmid from
    output_directory    Path to save fixed profiles
    remove_VGBs         Path to dataframe containing VGBs to drop
Description:
    This script finds all "_BAQLaVa_profile.txt" files in a directory and removes the set of plasmids found to be problematic in BAQLaVa manuscript (iHMP dataset)
"""

dropmgx = pd.read_csv(sys.argv[2], sep="\t")


def remove_plasmid_VGBs(loc1, loc2, plas):
    df1 = plas.copy()#.rename({'feature':'BAQLaVa VGB'})
    df1['featn'] = df1['feature']+"|nucleotide"
    df1['featt'] = df1['feature']+"|translated"
    df2 = pd.concat([df1.copy()[['feature']].rename(columns={'feature':'BAQLaVa VGB'}),
                     df1.copy()[['featn']].rename(columns={'featn':'BAQLaVa VGB'}),
                     df1.copy()[['featt']].rename(columns={'featt':'BAQLaVa VGB'})])

    for i in os.listdir(loc1):
        if "BAQLaVa_profile.txt" in i:
            if "BAQLaVa_profile_plasmid_removed.txt" not in i: 
                tempdf = pd.read_csv(loc1+i, sep="\t")
                tempdf2 = pd.merge(tempdf, df2.copy(), on='BAQLaVa VGB', how='left', indicator=True).query("_merge=='left_only'").drop(columns=['_merge'])
                tempdf2.to_csv(loc2 + "/" + i[:-4]+"_plasmid_removed.txt", sep="\t", index=False)

    return None

remove_plasmid_VGBs(sys.argv[1], sys.argv[3], dropmgx)
