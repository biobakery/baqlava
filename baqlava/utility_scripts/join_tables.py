# to run: python join_tables.py <DIRECTORY_WITH_BAQLAVA_TABLES> <DESIRED_PREFIX_TO_SAVE> 

import pandas as pd
import os
import sys

def main():
    retdf = pd.DataFrame({'BAQLaVa VGB':[]})
    for i in os.listdir(sys.argv[1]):
        if "_BAQLaVa_profile.txt" in i:
            name = i.split("_BAQLaVa_")[0] + "_Abundance-RPKs"
            tempdf = pd.read_csv(sys.argv[1] + i, sep="\t")
            tempdf = tempdf[~tempdf['BAQLaVa VGB'].str.contains("\|")]
            colname = tempdf.columns[1]
            tempdf = tempdf[['BAQLaVa VGB',colname]].rename(columns={colname:name})
            retdf = pd.merge(retdf, tempdf, on='BAQLaVa VGB', how='outer')
    a = retdf.fillna(0)
    a.to_csv(sys.argv[2] + "_BAQLaVa_VGB_table.tsv", sep="\t", index=False)
    return None

if __name__ == '__main__':
    main()
