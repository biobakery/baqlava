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
    default = "utility_files/idmap4.txt")

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

workflow.add_argument(
    name = "reconcile",
    desc = "script to reconcile reads",
    default = "reconcile_mapped_reads_v0.2.py")

args = workflow.parse_args()

############################
# function to run workflow #
############################

file_base = args.input.split("/")[-1].split(".")[0]
output_dir = args.output + "/" + file_base + "/"

workflow.add_task(
    "humann --input [depends[0]] --output [output_folder[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --protein-database [p_db] --threads [threads]",
    depends = [args.input],
    output_folder = [output_dir],
    targets = [output_dir + file_base + "_genefamilies.tsv"],
    n_db = args.nucdb,
    idx = args.nucindex,
    p_db = args.protdb,
    threads = args.threads)

workflow.add_task(
    "python [script] [depends[0]]" ,
    depends = [output_dir + file_base + "_genefamilies.tsv"],
    targets = [output_dir + file_base + "_baqlava_genefamilies.tsv"],
    output_folder = output_dir,
    script = args.reconcile)

workflow.go()
