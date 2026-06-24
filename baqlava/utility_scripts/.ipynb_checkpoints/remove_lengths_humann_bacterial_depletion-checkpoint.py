import sys
import os
import gzip
from pathlib import Path

#newname = os.path.split(sys.argv[1])[1].replace("_bowtie2_unaligned.fa", "_processed.fa")
#savefile = os.path.join(sys.argv[1], newname)

p = Path(sys.argv[1])
savefile = p.parent.parent / p.name.replace("_bowtie2_unaligned.fa","_processed.fa")


if sys.argv[1].endswith(".gz"):
    open_func=gzip.open
else:
    open_func=open

with open_func(sys.argv[1], "rt") as deplfasta:
    with open_func(savefile, "wt") as wfasta:
        for i in deplfasta:
            i = i.strip()
            if ">" in i:
                i = i.split("|")[:-1]
                i = "|".join(i)
                wfasta.write(i)
                wfasta.write("\n")
            else:
                wfasta.write(i)
                wfasta.write("\n")
