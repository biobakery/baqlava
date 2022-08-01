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
    desc = "nucleotide database folder to use",
    default = "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/nucleotide_database/")

workflow.add_argument(
    name = "nucindex",
    desc = "nucleotide annotation index to use",
    default = "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/additional_files/idmap3.txt")

workflow.add_argument(
    name = "input-extension",
    desc = "the input file extension",
    default = "fastq")

workflow.add_argument(
    name = "threads",
    desc="number of threads for knead_data to use",
    default=1)

workflow.add_argument(
    name = "protdb",
    desc = "protein database folder to use",
    default = "/n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/protein_database/")

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

humann_output_files_vir = args.output+"/humann_output_files_vir/"
temp_assembly_file = args.output+"/temp_assembly/"


workflow.add_task(
    "mkdir -p " + humann_output_files_vir,
    depends = fastq_files,
    targets = humann_output_files_vir,
    output_folder = args.output)

for base in file_bases:
    workflow.add_task(
        "humann --input [depends[0]] --output [output_folder[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --protein-database [p_db] --threads [threads]",
        depends = [file_path + base + "_cat.fastq"],
        output_folder = [humann_output_files_vir],
        targets = [humann_output_files_vir + base + "_cat_genefamilies.tsv"],
        n_db = args.nucdb,
        idx = args.nucindex,
        p_db = args.protdb,
        threads = args.threads)
    workflow.add_task(
        "python reconcile_mapped_reads.py [depends[0]]" ,
        depends = [humann_output_files_vir + base + "_cat_genefamilies.tsv" ],
        targets = [humann_output_files_vir + base + "_cat_baqlava_genefamilies.tsv"],
        output_folder = humann_output_files_vir)


workflow.go()


