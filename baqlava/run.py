"""
BAQLaVa: run module

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Viral Profiling"     #Update the description as needed
    )

### run.py [-i <folder_with_input_files> -o <output_directory>]

###############
# custom args #
###############

workflow.add_argument(
    name = "nucdb",
    desc = "nucleotide database to use",
    default = "####")

workflow.add_argument(
    name = "nucindex",
    desc = "nucleotide annotation index to use",
    default = "####")

workflow.add_argument(
    name = "input-extension",
    desc = "the input file extension",
    default = "fastq")

workflow.add_argument(
    name = "threads",
    desc="number of threads for knead_data to use",
    default=1)

args = workflow.parse_args()

############################
# function to run workflow #
############################

fastq_files = workflow.get_input_files(extension=args.input_extension)
file_path = args.input + "/"

def collect_and_cat_files(dir):
    bases = set()
    for file in dir:
        base = file.split("/")[-1]
        base = base.split("_")[0]
        bases.add(base)
    return list(bases)

file_bases = collect_and_cat_files(fastq_files)

for base in file_bases:
    workflow.add_task(
        "cat " + file_path + base + "_1.fastq " + file_path + base + "_2.fastq > " + file_path + base + "_cat.fastq",
        depends = fastq_files,
        targets = file_path + base + "_cat.fastq")

humann_output_files_bac = args.output+"/humann_output_files_bac/"
humann_output_files_vir = args.output+"/humann_output_files_vir/"
temp_assembly_file = args.output+"/temp_assembly/"

workflow.add_task(
    "mkdir -p " + humann_output_files_bac,
    depends = fastq_files,
    targets = humann_output_files_bac,
    output_folder = args.output)

workflow.add_task(
    "mkdir -p " + humann_output_files_vir,
    depends = fastq_files,
    targets = humann_output_files_vir,
    output_folder = args.output)

workflow.add_task(
    "mkdir -p " + temp_assembly_file,
    depends = fastq_files,
    targets = temp_assembly_file)

for base in file_bases:
    workflow.add_task(
        "humann --input [depends[0]] --output [output_folder[0]] --threads [threads]",
        depends = [file_path + base + "_cat.fastq", humann_output_files_bac],
        output_folder = [humann_output_files_bac],
        targets = [humann_output_files_bac + base +"_cat_pathabundance.tsv", humann_output_files_bac + base + "_cat_genefamilies.tsv"],
        threads = args.threads)
    workflow.add_task(
        "humann --input [depends[0]] --output [output_folder[0]] --bypass-nucleotide-index --nucleotide-database [db] --id-mapping [idx] --threads [threads]",
        depends = [file_path + base + "_cat.fastq", humann_output_files_bac],
        output_folder = [humann_output_files_vir],
        targets = [humann_output_files_vir + base + ".pathabundance.tsv"],
        db = args.nucdb,
        idx = args.nucindex,
        threads = args.threads)
    workflow.add_task(
        "mkdir -p [targets[0]]",
        depends = fastq_files,
        targets = [temp_assembly_file + base],
        output_folder = temp_assembly_file)
    workflow.add_task(
        "spades.py --meta -1 " + base + "_1.fastq -2 " + base + "_2.fastq -o [depends[0]]",
        depends = [temp_assembly_file + base],
        targets = [temp_assembly_file + base + "/contigs.fasta"],
        output_folder = temp_assembly_file + base)


workflow.go()



