import sys
import os
import gzip

# the bacterially depleted FASTA is always written gzipped to save space
newname = os.path.split( sys.argv[1] )[1].replace("_bowtie2_unaligned.fa", "_bacterial_depleted.fa.gz")
savefile = sys.argv[2] + newname

if sys.argv[1].endswith(".gz"):
    open_func=gzip.open
else:
    open_func=open

with open_func(sys.argv[1], "rt") as deplfasta:
    with gzip.open(savefile, "wt", compresslevel=6) as wfasta:
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
