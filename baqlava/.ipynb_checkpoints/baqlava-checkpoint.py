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
from pathlib import Path

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

class AbsolutePathConfigParser(configparser.ConfigParser):
    def get(self, section, option, **kwargs):
        # Get the raw value from the configuration file.
        value = super().get(section, option, **kwargs)
        # If the value is not absolute, prepend the library directory.
        if os.path.isabs(value):
            return value
        # If the value contains "//", assume it’s a URL and return it as is.
        if '//' in value:
            return value
        # If the value does not contain a slash at all, return it.
        if '/' not in value:
            return value
        # Otherwise, join the value with lib_dir to create an absolute path.
        return os.path.join(lib_dir, value)

lib_dir = os.path.dirname(os.path.abspath(__file__))
config_file=os.path.join(lib_dir,"configs/baqlava.cfg")
config = AbsolutePathConfigParser()
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
    default = config.getboolean('features','bypass_bacterial_depletion'))

workflow.add_argument(
    name = "bypass-nucleotide-search",
    desc = "Using this flag turns OFF nucleotide search. BAQLaVa will only use protein search in viral profiling.",
    action = 'store_true',
    default = config.getboolean("features", "bypass_nucleotide_search"))

workflow.add_argument(
    name = "bypass-translated-search",
    desc = "Using this flag turns OFF translated search. BAQLaVa will only use nucleotide search in viral profiling.",
    action = 'store_true',
    default = config.getboolean('features','bypass_translated_search'))

workflow.add_argument(
    name = "taxonomic-profile",
    desc = "Using this flag supplies a metaphlan generated bacterial species profile to baqlava to deplete bacterial reads rather than running the full metaphlan search.",
    default = config.getboolean('features','taxonomic_profile'))

workflow.add_argument(
    name = "proteome-length",
    desc = "Minimum length of proteome mapped to report in BAQLaVa output.",
    default = config.get('features','proteome_length'))

workflow.add_argument(
    name = "keep-tempfiles",
    desc = "Keep temporary mapping files. This is memory intensive.",
    action = 'store_true',
    default = config.get('features','keep_tempfiles'))

workflow.add_argument(
    name="genome-filtering",
    desc="Level of genome filtering to allow genomes to be represented in a VGB (genomes removed based on probability of plasmid). Options: no-filtering, default, conservative",
    default="default"
)

args = workflow.parse_args()

############ settings based on provided args:

GENOME_FILTERS = {
    "no-filtering": {
        "nucleotide": config.get('utility','idmap_nucl_NF'),
        "translated": config.get('utility','idmap_prot_NF'),
    },
    "default": {
        "nucleotide": config.get('utility','idmap_nucl_D'),
        "translated": config.get('utility','idmap_prot_D'),
    },
    "conservative": {
        "nucleotide": config.get('utility','idmap_nucl_C'),
        "translated": config.get('utility','idmap_prot_C'),
    },
}

nucleotide_idmapping = GENOME_FILTERS[args.genome_filtering]["nucleotide"]
translated_idmapping = GENOME_FILTERS[args.genome_filtering]["translated"]

############################# MAIN

def main():
    ############################
    # function to run workflow #
    ############################
    file_name = Path(args.input).name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fasta.gz", ".fa.gz", ".fasta", ".fa"]:
        if file_name.endswith(suffix):
            file_base = file_name.removesuffix(suffix)
            break
    else:
        file_base = Path(file_name).stem  # fallback for other extensions

    output_dir = Path(args.output)
    tempdir = output_dir / f"{file_base}_temp"
    baq_dir = output_dir / f"{file_base}_baqlava"

    logger.info(f"Output files will be written to: {output_dir}")
    message = f"Writing temp files to directory: {tempdir}"
    print("\n"+message+"\n")

    # create directories
    output_dir.mkdir(parents=True, exist_ok=True)
    tempdir.mkdir(parents=True, exist_ok=True)
    baq_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
     filename='log_file_name.log',
     level=logging.INFO,
     format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
     datefmt='%H:%M:%S')


    ########### File processing module ###########

    # !!!!!!!!add step to check if file has length appended to end

    if args.bypass_bacterial_depletion: # bypassing bacterial depletion
        # NO UPSTREAM BACTERIAL DEPLETION

            workflow.add_task(
                "scp [depends[0]] [targets[0]]",
                depends = [args.input],
                targets = [str(tempdir / f"{file_base}_processed.fa")],
                name = "Upstream file processing")

    elif not args.bypass_bacterial_depletion: # bacterial depletion
        # UPSTREAM BACTERIAL DEPLETION

        if not args.taxonomic_profile:
        # PROCEED WITHOUT METAPHLAN TAXONOMIC PROFILE
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-translated-search --threads [threads]",
                depends = [args.input],
                args = [tempdir],
                targets = [str(tempdir / f"{file_base}_humann_temp" / f"{file_base}_bowtie2_unaligned.fa")],
                threads = args.threads,
                name = "Running HUMAnN to deplete bacterial reads from file")

        else:
        # USE A METAPHLAN TAXONOMIC PROFILE TO AID BACTERIAL DEPLETION
            workflow.add_task(
                "humann --input [depends[0]] --output [args[0]] --bypass-translated-search --threads [threads] --taxonomic-profile [depends[1]]",
                depends = [args.input, args.taxonomic_profile],
                args = [tempdir],
                targets = [str(tempdir / f"{file_base}_humann_temp" / f"{file_base}_bowtie2_unaligned.fa")],
                threads = args.threads,
                name = "Running HUMAnN to deplete bacterial reads from file")

        workflow.add_task(
            "python [len_adj] [depends[0]]",
            depends = [str(tempdir / f"{file_base}_humann_temp" / f"{file_base}_bowtie2_unaligned.fa")],
            targets = [str(tempdir / f"{file_base}_processed.fa")],
            len_adj = os.path.abspath(args.lengthadjust),
            name = "Formatting bacterially depleted FASTA file")

        workflow.add_task(
            "scp [depends[0]] [targets[0]]",
            depends = [str(tempdir / f"{file_base}_processed.fa")],
            targets = [str(output_dir / f"{file_base}_bacterial_depleted.fa")],
            name = "Saving bacterially depleted FASTA file")

    ########### BAQLaVa NUCLEOIDE & TRANSLATED SEARCH ###########

    if not args.bypass_nucleotide_search: # bypass nucleotide != True

        # FIRST RUN & CALCULATE AT 25% COVERAGE:
        workflow.add_task(
            "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --count-normalization 'Adjusted RPKs' --nucleotide-subject-coverage-threshold 25",
            depends = [str(tempdir / f"{file_base}_processed.fa")],
            args = [baq_dir, str(f"{file_base}_nucleotide1")],
            targets = [str(baq_dir / f"{file_base}_nucleotide1_2_genefamilies.tsv")],
            n_db = os.path.abspath(args.nucdb),
            idx = os.path.abspath(nucleotide_idmapping),
            threads = args.threads,
            name = "Running BAQLaVa Nucleotide Search")

        workflow.add_task(
            "scp [depends[0]] [targets[0]]",
            depends = [str(baq_dir / f"{file_base}_nucleotide1_2_genefamilies.tsv")],
            targets = [str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv")],
            name = "File Processing")

        # USE FIRST RUN TO CALCULATE AT 50% COVERAGE:
        workflow.add_task(
            "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --count-normalization 'Adjusted RPKs' --nucleotide-subject-coverage-threshold 50 --resume",
            depends = [str(tempdir / f"{file_base}_processed.fa"), str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv")],
            args = [baq_dir, str(f"{file_base}_nucleotide2")],
            targets = [str(baq_dir / f"{file_base}_nucleotide2_2_genefamilies.tsv")],
            n_db = os.path.abspath(args.nucdb),
            idx = os.path.abspath(nucleotide_idmapping),
            threads = args.threads,
            name = "Calculating Marker Coverage & Abundance")

        workflow.add_task(
            "scp [depends[0]] [targets[0]]",
            depends = [str(baq_dir / f"{file_base}_nucleotide2_2_genefamilies.tsv"), str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv")],
            targets = [str(baq_dir / f"{file_base}_nucleotide_50_genefamilies.tsv")],
            name = "File Processing")

    if not args.bypass_translated_search: # bypass translated != True

        workflow.add_task(
            "humann --input [depends[0]] --output [args[0]] --id-mapping [idx] --protein-database [p_db] --threads [threads] --bypass-nucleotide-search --output-basename [args[1]] --count-normalization 'Adjusted RPKs' --translated-subject-coverage-threshold 50",
            depends = [args.input],
            args = [baq_dir, str(f"{file_base}_translated")],
            targets = [str(baq_dir / f"{file_base}_translated_2_genefamilies.tsv")],
            p_db = os.path.abspath(args.protdb),
            idx = os.path.abspath(translated_idmapping),
            threads = args.threads,
            name = "Running BAQLaVa Translated Search")


    ########### RECONCILING GENEFAMILIES.TSV FILES BASED ON PROCESSES RUN ###########

    if not args.bypass_nucleotide_search and args.bypass_translated_search: # nucleotide only, no translated
  # if args.bypass_nucleotide_search == 'False' and args.bypass_translated_search == True:

        workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [args[0]] [depends[2]] [proteomelen] [targets[0]] [targets[1]] [args[0]]",
            script = args.reconcile_mapped_script,
            mode = "1",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv"), str(baq_dir / f"{file_base}_nucleotide_50_genefamilies.tsv"), args.input],
            targets = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt"), str(output_dir / f"{file_base}_tempfile_markers.txt")],
            name = "Making BAQLaVa Viral Profile")

    elif args.bypass_nucleotide_search and not args.bypass_translated_search: # translated only, no nucleotide
#        elif args.bypass_nucleotide_search == True and args.bypass_translated_search == 'False':

        workflow.add_task("python [script] [mode] [args[0]] [args[0]] [depends[0]] [args[0]] [proteomelen] [targets[0]] [args[0]] [targets[1]]",
            script = args.reconcile_mapped_script,
            mode = "2",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [str(baq_dir / f"{file_base}_translated_2_genefamilies.tsv")],
            targets = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt"), str(output_dir / f"{file_base}_tempfile_proteins.txt")],
            name = "Making BAQLaVa Viral Profile")

    elif not args.bypass_nucleotide_search and not args.bypass_translated_search: # nucleotide + translated (full workflow)
#        elif args.bypass_nucleotide_search == 'False' and args.bypass_translated_search == 'False':

            workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [depends[2]] [depends[3]] [proteomelen] [targets[0]] [targets[1]] [targets[2]]",
            script = args.reconcile_mapped_script,
            mode = "3",
            proteomelen = args.proteome_length,
            depends = [str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv"), str(baq_dir / f"{file_base}_nucleotide_50_genefamilies.tsv"), str(baq_dir / f"{file_base}_translated_2_genefamilies.tsv"), args.input],
            targets = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt"), str(output_dir / f"{file_base}_tempfile_markers.txt"), str(output_dir / f"{file_base}_tempfile_proteins.txt")],
            name = "Making BAQLaVa Viral Profile")

    else:
        print("ERROR: Incorrect reconciliation task prompted")

#    elif args.bypass_bacterial_depletion == 'False':
#    # UPSTREAM BACTERIAL DEPLETION

#        if args.taxonomic_profile == 'False':
#        # PROCEED WITHOUT METAPHLAN TAXONOMIC PROFILE

#            workflow.add_task(
#                "humann --input [depends[0]] --output [args[0]] --bypass-translated-search --threads [threads]",
#                depends = [args.input],
#                args = [tempdir],
#                targets = [tempdir + "/" + file_base + "_humann_temp/" + file_base + "_bowtie2_unaligned.fa"],
#                threads = args.threads,
#                name = "Running HUMAnN to depete bacterial reads from file")

#        else:
#        # USE A METAPHLAN TAXONOMIC PROFILE TO AID BACTERIAL DEPLETION

#            workflow.add_task(
#                "humann --input [depends[0]] --output [args[0]] --bypass-translated-search --threads [threads] --taxonomic-profile [depends[1]]",
#                depends = [args.input, args.taxonomic_profile],
#                args = [tempdir],
#                targets = [tempdir + "/" + file_base + "_humann_temp/" + file_base + "_bowtie2_unaligned.fa"],
#                threads = args.threads,
#                name = "Running HUMAnN to depete bacterial reads from file")

#        workflow.add_task(
#            "python [len_adj] [depends[0]] [args[0]]",
#            depends = [tempdir + "/" + file_base + "_humann_temp/" + file_base + "_bowtie2_unaligned.fa"],
#            args = [output_dir],
#            targets = [output_dir + file_base + "_bacterial_depleted.fa"],
#            len_adj = os.path.abspath(args.lengthadjust),
#            name = "Formatting bacterially depleted FASTA file")

#        # BAQLAVA VIRAL PROFILING:
#        if args.bypass_nucleotide_search == True:
#            # BYPASS NUCLEOTIDE SEARCH:
#            pass
#        else:

#            # NUCLEOTIDE SEARCH:
#            # FIRST RUN & CALCULATE AT 25% COVERAGE:
#            workflow.add_task(
#                "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --nucleotide-subject-coverage-threshold 25",
#                depends = [output_dir + file_base + "_bacterial_depleted.fa"],
#                args = [baq_dir, file_base + "_bacterial_depleted_nucleotide"],
#                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv"],
#                n_db = os.path.abspath(args.nucdb),
#                idx = os.path.abspath(nucleotide_idmapping),
#                threads = args.threads,
#                name = "Running BAQLaVa Nucleotide Search")

#            workflow.add_task(
#                "mv [depends[0]] [targets[0]]",
#                depends = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv"],
#                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv"])

#            # USE FIRST RUN TO CALCULATE AT 50% COVERAGE:
#            workflow.add_task(
#                "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --nucleotide-subject-coverage-threshold 50 --resume",
#                depends = [output_dir + file_base + "_bacterial_depleted.fa", baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv"],
#                args = [baq_dir, file_base + "_bacterial_depleted_nucleotide"],
#                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv"],
#                n_db = os.path.abspath(args.nucdb),
#                idx = os.path.abspath(nucleotide_idmapping),
#                threads = args.threads,
#                name = "Calculating Marker Coverage & Abundance")

#            workflow.add_task(
#                "mv [depends[0]] [targets[0]]",
#                depends = [baq_dir + file_base + "_bacterial_depleted_nucleotide_genefamilies.tsv", baq_dir + file_base + "_bacterial_depleted_nucleotide_25_genefamilies.tsv"],
#                targets = [baq_dir + file_base + "_bacterial_depleted_nucleotide_50_genefamilies.tsv"])

#        if args.bypass_translate

#        else:
#            print("ERROR: Incorrect reconciliation task prompted")

#    else:
#        print('ERROR')


#    if args.keep_tempfiles == True:
#        pass
#    else:
#        workflow.add_task(
#        "rm -r [args[0]]",
#        depends = [output_dir + file_base + "_BAQLaVa_profile.txt"],
#        args = [baq_dir])

#    workflow.add_task(
#    "rm -r [args[0]]",
#    depends = [output_dir + file_base + "_BAQLaVa_profile.txt"],
#    args = [tempdir])

    workflow.go()


if __name__ == '__main__':
    main()
