#!/usr/bin/env python
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
import sys
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable
import tempfile
import logging

# name global logging instance
logger=logging.getLogger(__name__)

# Config Parsers
# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser

try:
    from baqlava import utility_scripts
except ImportError:
    sys.exit("ERROR: Unable to find the baqlava python package." +
        " Please check your install.")

# Config Parsers
# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser
try:
    import ConfigParser as configparser
except ImportError:
    import configparser

config = configparser.ConfigParser()
install_folder=os.path.dirname(os.path.realpath(__file__))
config_file=os.path.join(install_folder,"configs/baqlava.cfg")
config.read(config_file)

# Setting the version of the workflow and short description

workflow = Workflow(
    version=config.get('metadata','version'),                    #Update the version as needed
    description="Viral Profiling"     #Update the description as needed
    )

###############
# custom args #
###############
workflow.add_argument(
    name = "nucdb",
    desc = "nucleotide database folder to use",
    default = config.get('database','nucleotide_db'))

workflow.add_argument(
    name = "nucindex",
    desc = "nucleotide annotation index to use",
    default = config.get('utility','idmap_nucl'))

workflow.add_argument(
    name = "protindex",
    desc = "protein annotation index to use",
    default = config.get('utility','idmap_prot'))

workflow.add_argument(
    name = "threads",
    desc="number of threads for each process to use",
    default=config.get('computation','threads'))

workflow.add_argument(
    name = "protdb",
    desc = "protein database folder to use",
    default = config.get('database','translated_db'))

workflow.add_argument(
    name = "lengthadjust",
    desc = "location of python script used to adjust HUMAnN temp file",
    default = config.get('features','humann_length_adjust'))

workflow.add_argument(
    name = "reconcile-mapped-script",
    desc = "location of python script used to run final reconciliation",
    default = config.get('code','reconcile_mapped_script'))

workflow.add_argument(
    name = "bypass-bacterial-depletion",
    desc = "Using this flag turns OFF bactieral depletion. Input file will be immediately profiled by BAQLaVa.",
    action = 'store_true',
    default = config.get('features','bypass_bacterial_depletion'))

workflow.add_argument(
    name = "bypass-nucleotide-search",
    desc = "Using this flag turns OFF nucleotide search. BAQLaVa will only use protein search in viral profiling.",
    action = 'store_true',
    default = config.get('features','bypass_nucleotide_search'))

workflow.add_argument(
    name = "bypass-translated-search",
    desc = "Using this flag turns OFF translated search. BAQLaVa will only use nucleotide search in viral profiling.",
    action = 'store_true',
    default = config.get('features','bypass_translated_search'))

workflow.add_argument(
    name = "taxonomic-profile",
    desc = "Using this flag supplies a metaphlan generated bacterial species profile to baqlava to deplete bacterial reads rather than running the full metaphlan search.",
    default = config.get('features','taxonomic_profile'))

workflow.add_argument(
    name = "proteome-length",
    desc = "Minimum length of proteome mapped to report in BAQLaVa output.",
    default = config.get('features','proteome_length'))

workflow.add_argument(
    name = "keep-tempfiles",
    desc = "Keep temporary mapping files. This is memory intensive.",
    action = 'store_true',
    default = config.get('features','keep_tempfiles'))

args = workflow.parse_args()


def main():
    ############################
    # function to run workflow #
    ############################
    file_base = args.input.split("/")[-1].split(".")[0]
    output_dir = args.output + "/"
    tempdir = args.output + "/" + file_base + "_temp/"
    baq_dir = args.output + "/" + file_base + "_baqlava/"

    logger.info("Output files will be written to: " + output_dir)
    message="Writing temp files to directory: " + tempdir
    print("\n"+message+"\n")

    # create directories
    os.system("mkdir " + args.output)
    os.system("mkdir " + tempdir)
    os.system("mkdir " + baq_dir)

    logging.basicConfig(
     filename='log_file_name.log',
     level=logging.INFO,
     format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
     datefmt='%H:%M:%S')

    # set the name of the log file
    #log_file=os.path.join(output_dir,config.file_basename+"_0.log")

    # change file name if set
    #if args.o_log:
    #    log_file=args.o_log

    # configure the logger
    #logging.basicConfig(filename=log_file,format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
    #    level=getattr(logging,args.log_level), filemode='w', datefmt='%m/%d/%Y %I:%M:%S %p')


    if args.bypass_bacterial_depletion == True:
        # NO UPSTREAM BACTERIAL DEPLETION

        if args.bypass_nucleotide_search == True:
            # BYPASSING NUCLEOTIDE SEARCH
            pass
        else:

            # NUCLEOTIDE SEARCH:
            # FIRST RUN & CALCULATE AT 25% COVERAGE:
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --nucleotide-subject-coverage-threshold 25",
                depends = [args.input],
                args = [baq_dir, file_base + "_nucleotide"],
                targets = [baq_dir + file_base + "_nucleotide_genefamilies.tsv"],
                n_db = os.path.abspath(args.nucdb),
                idx = os.path.abspath(args.nucindex),
                threads = args.threads,
                name = "Running BAQLaVa Nucleotide Search")

            workflow.add_task(
                "mv [depends[0]] [targets[0]]",
                depends = [baq_dir + file_base + "_nucleotide_genefamilies.tsv"],
                targets = [baq_dir + file_base + "_nucleotide_25_genefamilies.tsv"])

            # USE FIRST RUN TO CALCULATE AT 50% COVERAGE:
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --nucleotide-subject-coverage-threshold 50 --resume",
                depends = [args.input, baq_dir + file_base + "_nucleotide_25_genefamilies.tsv"],
                args = [baq_dir, file_base + "_nucleotide"],
                targets = [baq_dir + file_base + "_nucleotide_genefamilies.tsv"],
                n_db = os.path.abspath(args.nucdb),
                idx = os.path.abspath(args.nucindex),
                threads = args.threads,
                name = "Calculating Marker Coverage & Abundance")

            workflow.add_task(
                "mv [depends[0]] [targets[0]]",
                depends = [baq_dir + file_base + "_nucleotide_genefamilies.tsv", baq_dir + file_base + "_nucleotide_25_genefamilies.tsv"],
                targets = [baq_dir + file_base + "_nucleotide_50_genefamilies.tsv"])

        if args.bypass_translated_search == True:
       	    # BYPASSING TRANSLATED SEARCH
            pass
       	else:
            # TRANSLATED SEARCH:
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --id-mapping [idx] --protein-database [p_db] --threads [threads] --bypass-nucleotide-search --output-basename [args[1]] --translated-subject-coverage-threshold 50",
                depends = [args.input],
                args = [baq_dir, file_base + "_translated"],
                targets = [baq_dir + file_base + "_translated_genefamilies.tsv"],
                p_db = os.path.abspath(args.protdb),
                idx = os.path.abspath(args.protindex),
                threads = args.threads,
                name = "Running BAQLaVa Translated Search")

        # NOW RECONCILE GENEFAMILIES.TSV FILES BASED ON PROCESESS RUN:
        if args.bypass_nucleotide_search == 'False' and args.bypass_translated_search == True:

            workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [args[0]] [depends[2]] [proteomelen] [targets[0]] [targets[1]] [args[0]]",
            script = args.reconcile_mapped_script,
            mode = "1",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [baq_dir + file_base + "_nucleotide_25_genefamilies.tsv", baq_dir + file_base + "_nucleotide_50_genefamilies.tsv", args.input],
            targets = [output_dir + file_base + "_BAQLaVa_profile.txt", output_dir + file_base + "_tempfile_markers.txt"],
            name = "Making BAQLaVa Viral Profile")

        elif args.bypass_nucleotide_search == True and args.bypass_translated_search == 'False':

            workflow.add_task("python [script] [mode] [args[0]] [args[0]] [depends[0]] [args[0]] [proteomelen] [targets[0]] [args[0]] [targets[1]]",
            script = args.reconcile_mapped_script,
            mode = "2",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [baq_dir + file_base + "_translated_genefamilies.tsv"],
            targets = [output_dir + file_base + "_BAQLaVa_profile.txt", output_dir + file_base + "_tempfile_proteins.txt"],
            name = "Making BAQLaVa Viral Profile")

        elif args.bypass_nucleotide_search == 'False' and args.bypass_translated_search == 'False':

            workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [depends[2]] [depends[3]] [proteomelen] [targets[0]] [targets[1]] [targets[2]]",
            script = args.reconcile_mapped_script,
            mode = "3",
            proteomelen = args.proteome_length,
            depends = [baq_dir + file_base + "_nucleotide_25_genefamilies.tsv", baq_dir + file_base + "_nucleotide_50_genefamilies.tsv", baq_dir + file_base + "_translated_genefamilies.tsv", args.input],
            targets = [output_dir + file_base + "_BAQLaVa_profile.txt", output_dir + file_base + "_tempfile_markers.txt", output_dir + file_base + "_tempfile_proteins.txt"],
            name = "Making BAQLaVa Viral Profile")

        else:
            print("ERROR: Incorrect reconciliation task prompted")

    elif args.bypass_bacterial_depletion == 'False':
    # UPSTREAM BACTERIAL DEPLETION

        if args.taxonomic_profile == 'False':
        # PROCEED WITHOUT METAPHLAN TAXONOMIC PROFILE

            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-translated-search --threads [threads]",
                depends = [args.input],
                args = [tempdir],
                targets = [tempdir + "/" + file_base + "_humann_temp/" + file_base + "_bowtie2_unaligned.fa"],
                threads = args.threads,
                name = "Running HUMAnN to depete bacterial reads from file")

        else:
        # USE A METAPHLAN TAXONOMIC PROFILE TO AID BACTERIAL DEPLETION

            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-translated-search --threads [threads] --taxonomic-profile [depends[1]]",
                depends = [args.input, args.taxonomic_profile],
                args = [tempdir],
                targets = [tempdir + "/" + file_base + "_humann_temp/" + file_base + "_bowtie2_unaligned.fa"],
                threads = args.threads,
                name = "Running HUMAnN to depete bacterial reads from file")

        workflow.add_task(
            "python [len_adj] [depends[0]] [args[0]]",
            depends = [tempdir + "/" + file_base + "_humann_temp/" + file_base + "_bowtie2_unaligned.fa"],
            args = [output_dir],
            targets = [output_dir + file_base + "_bacterial_depleted.fa"],
            len_adj = os.path.abspath(args.lengthadjust),
            name = "Formatting bacterially depleted FASTA file")

        # BAQLAVA VIRAL PROFILING:
        if args.bypass_nucleotide_search == True:
            # BYPASS NUCLEOTIDE SEARCH:
            pass
        else:

            # NUCLEOTIDE SEARCH:
            # FIRST RUN & CALCULATE AT 25% COVERAGE:
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --nucleotide-subject-coverage-threshold 25",
                depends = [output_dir + file_base + "_bacterial_depleted.fa"],
                args = [baq_dir, file_base + "_bacterial_depleted_nucleotide"],
                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv"],
                n_db = os.path.abspath(args.nucdb),
                idx = os.path.abspath(args.nucindex),
                threads = args.threads,
                name = "Running BAQLaVa Nucleotide Search")

            workflow.add_task(
                "mv [depends[0]] [targets[0]]",
                depends = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv"],
                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv"])

            # USE FIRST RUN TO CALCULATE AT 50% COVERAGE:
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --nucleotide-subject-coverage-threshold 50 --resume",
                depends = [output_dir + file_base + "_bacterial_depleted.fa", baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv"],
                args = [baq_dir, file_base + "_bacterial_depleted_nucleotide"],
                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv"],
                n_db = os.path.abspath(args.nucdb),
                idx = os.path.abspath(args.nucindex),
                threads = args.threads,
                name = "Calculating Marker Coverage & Abundance")

            workflow.add_task(
                "mv [depends[0]] [targets[0]]",
                depends = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv", baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv"],
                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_50_genefamilies.tsv"])

        if args.bypass_translated_search == True:
            # BYPASS TRANSLATED SEARCH
            pass
        else:
            # TRANSLATED SEARCH:
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --id-mapping [idx] --protein-database [p_db] --threads [threads] --bypass-nucleotide-search --output-basename [args[1]] --translated-subject-coverage-threshold 50",
                depends = [output_dir + file_base + "_bacterial_depleted.fa"],
                args = [baq_dir, file_base + "_bacterial_depleted_translated"],
                targets = [baq_dir + file_base + "_bacterial_depleted_translated_genefamilies.tsv"],
                p_db = os.path.abspath(args.protdb),
                idx = os.path.abspath(args.protindex),
                threads = args.threads,
                name = "Running BAQLaVa Translated Search")

        # NOW RECONCILE GENEFAMILIES.TSV FILES BASED ON PROCESESS RUN:

        if args.bypass_nucleotide_search == 'False' and args.bypass_translated_search == True:

            workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [args[0]] [depends[2]] [proteomelen] [targets[0]] [targets[1]] [args[0]]",
            script = args.reconcile_mapped_script,
            mode = "1",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv", baq_dir + file_base + "_bacterial_depleted_nucleotide_50_genefamilies.tsv", args.input],
            targets = [output_dir + file_base + "_BAQLaVa_profile.txt", output_dir + file_base + "_tempfile_markers.txt"],
            name = "Making BAQLaVa Viral Profile")

        elif args.bypass_nucleotide_search == True and args.bypass_translated_search == 'False':

            workflow.add_task("python [script] [mode] [args[0]] [args[0]] [depends[0]] [args[0]] [proteomelen] [targets[0]] [args[0]] [targets[1]]",
            script = args.reconcile_mapped_script,
            mode = "2",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [baq_dir + file_base + "_bacterial_depleted_translated_genefamilies.tsv"],
            targets = [output_dir + file_base + "_BAQLaVa_profile.txt", output_dir + file_base + "_tempfile_proteins.txt"],
            name = "Making BAQLaVa Viral Profile")

        elif args.bypass_nucleotide_search == 'False' and args.bypass_translated_search == 'False':

            workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [depends[2]] [depends[3]] [proteomelen] [targets[0]] [targets[1]] [targets[2]]",
            script = args.reconcile_mapped_script,
            mode = "3",
            proteomelen = args.proteome_length,
            depends = [baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv", baq_dir + file_base + "_bacterial_depleted_nucleotide_50_genefamilies.tsv", baq_dir + file_base + "_bacterial_depleted_translated_genefamilies.tsv", args.input],
            targets = [output_dir + file_base + "_BAQLaVa_profile.txt", output_dir + file_base + "_tempfile_markers.txt", output_dir + file_base + "_tempfile_proteins.txt"],
            name = "Making BAQLaVa Viral Profile")


        else:
            print("ERROR: Incorrect reconciliation task prompted")

    else:
        print('ERROR')


    workflow.add_task(
    "rm -r [args[0]]",
    depends = [output_dir + file_base + "_BAQLaVa_profile.txt"],
    args = [baq_dir])

    if args.keep_tempfiles == True:
        pass
    else:
        workflow.add_task(
        "rm -r [args[0]]",
        depends = [output_dir + file_base + "_BAQLaVa_profile.txt"],
        args = [tempdir])

    workflow.go()


if __name__ == '__main__':
    main()
