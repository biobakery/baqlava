import sys
import os
import gzip

savefile = sys.argv[1].replace(".fa", "_lengthremoved.fa")

if sys.argv[1].endswith(".gz"):
    open_func=gzip.open
else:
    open_func=open

with open_func(sys.argv[1], "rt") as deplfasta:
    with open_func(savefile, "wt") as wfasta:
        for i in deplfasta:
            i = i.strip()
            if ">" in i:
                #print(i)
                i = i.split("|")[:-1]
                i = "|".join(i)
                wfasta.write(i)
                wfasta.write("\n")
            else:
                wfasta.write(i)
                wfasta.write("\n")
