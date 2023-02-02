import sys
import os

savefile = sys.argv[1].replace(".fa", ".lengthremoved.fa")

with open(sys.argv[1], "r") as deplfasta:
    with open(savefile, "w") as wfasta:
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
