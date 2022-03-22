import pandas as pd
import sys

def collect_reads(file):
    reads = []
    with open(file, "r") as fasta:
        for line in fasta:
            line=line.strip()
            if line[0] == ">":
                reads.append(line)
    return reads

def find_unmapped_intersection(lis_H, lis_V):
    HUMAnN = set(lis_H)
    viral = set(lis_V)
    unmapped_reads = HUMAnN.intersection(viral)
    return unmapped_reads

def get_saveout_name(sysargv1):
    name = str(sysargv1)
    name = name.split("/")[-1]
    name = name.replace('_diamond_unaligned.fa', '')
    name = name + '_combined_unaligned.fa'
    return name

def write_out_unaligned(file, unaligned_reads, name):
    keeper_counter = 0
    with open(file, "r") as fasta:
        with open(name, "w") as writeout:
            for line in fasta:
                line = line.strip()
                if line[0] == ">":
                    if line in unaligned_reads:
                        writeout.write(line)
                        writeout.write("\n")
                        keeper_counter = 1
                    else:
                        keeper_counter = 0
                else:
                    if keeper_counter == 1:
                        writeout.write(line)
                        writeout.write("\n")
                        keeper_counter = 0
                    else:
                        pass
    return None


HUMAnN_unmapped = collect_reads(sys.argv[1])
viral_unmapped = collect_reads(sys.argv[1])

unmapped = find_unmapped_intersection(HUMAnN_unmapped,viral_unmapped)

name = get_saveout_name(sys.argv[1])

write_out_unaligned(sys.argv[1], unmapped, name)
