import pandas as pd
import os
import sys

nuc_tax_file = pd.read_csv(sys.argv[1], sep="\t", index_col='Unnamed: 0', low_memory=False)

def species_formatting(df):
    df1 = df.copy()
    df1['segment_group'] = df1['segment_group'].fillna(df1['VGB'])
    df2 = df1[['Species','segment_group']].drop_duplicates()
    format_names = []
    for i in range(len(df2)):
        if 'segment' in df2.iloc[i]['segment_group']:
            name = df2.iloc[i]['Species']
            name2 = "_".join(name.split(" ")).upper()
            format_names.append(name2 + ".SEGMENTED_GROUP_" + df2.iloc[i]['segment_group'].split("_")[2])
        elif 'VGB' in df2.iloc[i]['segment_group']:
            VGB = df2.iloc[i]['Species']
            if len(VGB.split("_")) == 2:
                name = "UNCLASSIFIED_SPECIES."+df2.iloc[i]['Species']
                format_names.append(name)
            else:
                if 'VGB' in VGB:
                    vgbnumber = "_".join(VGB.split("_")[0:2])
                    speciesgroup = "_".join(VGB.split("_")[2:]).upper()
                    name = speciesgroup + "." + vgbnumber
                    format_names.append(name)
                else:
                    name = "_".join(VGB.split(" ")).upper() + "." + df2.iloc[i]['segment_group']
                    format_names.append(name)
    df2['Species_formatted'] = format_names
    return df2.rename(columns={'segment_group':'VGB'})

species_formatting(nuc_tax_file).to_csv(sys.argv[2], sep="\t")
