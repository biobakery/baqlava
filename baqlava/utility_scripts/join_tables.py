import pandas as pd
import os
import sys

def main():
    retdf = pd.DataFrame({'Species':[]})
    for i in os.listdir(sys.argv[1]):
        if "_BAQLaVa_profile.tsv" in i:
            name = i.split("_BAQLaVa_")[0] + "_Abundance-RPKs"
            colname = i.split("_BAQLaVa_")[0] + "_bowtie2_unaligned_lengthremoved_Abundance-RPKs"
            tempdf = pd.read_csv(sys.argv[1] + i, sep="\t").query("Database=='TOTAL'").drop(columns=['Database']).rename(columns={colname:name})
            retdf = pd.merge(retdf, tempdf, on='Species', how='outer')
    a = retdf.fillna(0)
    a.to_csv(sys.argv[2] + "BAQLaVa_VGB_table.tsv", sep="\t", index=False)
    return None

if __name__ == '__main__':
    main()
