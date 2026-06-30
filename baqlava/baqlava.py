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
import re
import subprocess
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
    name = "bacterial-depletion-script",
    desc = "location of python script used to run bacterial depletion (fails gracefully and falls back to the original input)",
    default = config.get('code','bacterial_depletion_script'))

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
    default = config.getboolean('features','keep_tempfiles'))

workflow.add_argument(
    name="genome-filtering",
    desc="Level of genome filtering to allow genomes to be represented in a VGB (genomes removed based on probability of plasmid). Options: no-filtering, default, conservative",
    default="default")

workflow.add_argument(
    name = "humann-depleted-fasta",
    desc = "This flag is designed to be used within the biobakery workflows. However, this flag can additionally be used if processing samples similar to the biobakery workflow - that is, if providing BAQLaVa with a fasta file that has already been processed by HUMAnN to deplete bacterial reads (HUMAnN output file named <SAMPLE>_bowtie2_unaligned.fa). If using this flag, --bypass-bacterial-depletion does NOT need to be additionally provided.",
    action = 'store_true',
    default = config.getboolean('features','humann_depleted_fasta'))

workflow.add_argument(
    name = "humann-passthrough-parameters-nucleotide",
    desc = "Provide nucleotide search parameters to HUMAnN - to maintain parameter choice as validated in BAQLaVa even if HUMAnN parameters change.",
    default = config.get('features','humann_passthrough_parameters_nucleotide'))

workflow.add_argument(
    name = "humann-passthrough-parameters-translated",
    desc = "Provide translated search parameters to HUMAnN - to maintain parameter choice as validated in BAQLaVa even if HUMAnN parameters change.",
    default = config.get('features','humann_passthrough_parameters_translated'))

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

############################# HUMAnN COMPATIBILITY CHECK

# BAQLaVa v1.2.0 requires HUMAnN 4 (H4). H4 uses the standard MetaPhlAn (MPA)
# database for its bacterial-depletion prescreen and performs its own MetaPhlAn
# DB version check at runtime, so BAQLaVa only needs to confirm that H4 (not an
# older HUMAnN) is the version on the PATH.
MINIMUM_HUMANN_MAJOR_VERSION = 4

def check_humann_version(minimum_major=MINIMUM_HUMANN_MAJOR_VERSION):
    """Verify HUMAnN >= 4.0 (H4) is installed; exit with an error otherwise."""
    try:
        proc = subprocess.run(
            ["humann", "--version"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, check=True)
    except FileNotFoundError:
        sys.exit("ERROR: HUMAnN was not found on your PATH. BAQLaVa requires "
                 "HUMAnN >= 4.0 (H4). Please install HUMAnN 4 and try again.")
    except subprocess.CalledProcessError as e:
        sys.exit("ERROR: Unable to determine the HUMAnN version "
                 "('humann --version' failed). BAQLaVa requires HUMAnN >= 4.0 "
                 "(H4).\n" + (e.stderr or ""))

    version_output = ((proc.stdout or "") + (proc.stderr or "")).strip()

    # e.g. "humann v4.0.0.alpha.2" -> major version 4
    match = re.search(r'humann\s+v?(\d+)', version_output, re.IGNORECASE)
    if not match:
        match = re.search(r'(\d+)', version_output)
    if not match:
        sys.exit("ERROR: Unable to parse the HUMAnN version from "
                 "'humann --version' output: '" + version_output + "'. "
                 "BAQLaVa requires HUMAnN >= 4.0 (H4).")

    major = int(match.group(1))
    if major < minimum_major:
        sys.exit("ERROR: BAQLaVa requires HUMAnN >= 4.0 (H4), but found "
                 "'" + version_output + "'. HUMAnN 4 uses the standard MetaPhlAn "
                 "(MPA) database required for BAQLaVa's bacterial-depletion step. "
                 "Please upgrade to HUMAnN 4.")

    logger.info("Verified HUMAnN version (H4+): " + version_output)
    print("\nVerified HUMAnN version (H4+): " + version_output + "\n")
    return version_output

############################# VERSION & DATABASE REPORTING

# Common MetaPhlAn DB tag style, e.g. vOct22_CHOCOPhlAnSGB_202403
MPA_DB_TAG = re.compile(r'v[A-Z][a-z]{2}\d{2}_[A-Za-z0-9]+_\d{6}')

def _run_version_cmd(cmd):
    """Return combined stdout+stderr of a command, or None if it cannot run."""
    try:
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)
    except (FileNotFoundError, OSError):
        return None
    if proc.returncode != 0:
        return None
    output = ((proc.stdout or "") + (proc.stderr or "")).strip()
    return output if output else None

def log_versions(humann_version=None):
    """Record BAQLaVa, sub-tool (HUMAnN, MetaPhlAn), MPA database, HUMAnN
    database, and BAQLaVa database versions to the logfile and stdout."""
    lines = ["BAQLaVa run environment - versions & databases:"]

    # BAQLaVa
    lines.append("  BAQLaVa version: " + config.get('metadata', 'version'))

    # HUMAnN (reuse the version captured at startup when available)
    if humann_version is None:
        humann_version = _run_version_cmd(["humann", "--version"])
    lines.append("  HUMAnN version: " + (humann_version or "unable to determine"))

    # MetaPhlAn tool version + MPA (MetaPhlAn) database version
    mpa_output = _run_version_cmd(["metaphlan", "--version"])
    if mpa_output:
        lines.append("  MetaPhlAn version: " + mpa_output.replace("\n", " ").strip())
        db_match = MPA_DB_TAG.search(mpa_output)
        lines.append("  MPA database: " + (
            db_match.group(0) if db_match
            else "unable to determine from 'metaphlan --version'"))
    else:
        lines.append("  MetaPhlAn version: unable to determine")
        lines.append("  MPA database: unable to determine")

    # HUMAnN databases (nucleotide ChocoPhlAn, protein UniRef, utility mapping)
    humann_cfg = _run_version_cmd(["humann_config", "--print"])
    db_folder_lines = []
    if humann_cfg:
        db_folder_lines = [ln.strip() for ln in humann_cfg.splitlines()
                           if "database_folders" in ln]
    if db_folder_lines:
        lines.append("  HUMAnN databases:")
        for ln in db_folder_lines:
            lines.append("    " + ln)
    else:
        lines.append("  HUMAnN databases: unable to determine")

    # BAQLaVa databases
    lines.append("  BAQLaVa nucleotide database: " + os.path.abspath(args.nucdb))
    lines.append("  BAQLaVa protein database: " + os.path.abspath(args.protdb))
    lines.append("  BAQLaVa genome-filtering level: " + str(args.genome_filtering))

    report = "\n".join(lines)
    logger.info(report)
    print("\n" + report + "\n")

############################# MAIN

def main():
    ############################
    # function to run workflow #
    ############################

    # Confirm HUMAnN 4 (H4) is available before doing any work.
    humann_version = check_humann_version()

    file_name = Path(args.input).name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fasta.gz", ".fa.gz", ".fasta", ".fa"]:
        if file_name.endswith(suffix):
            file_base = file_name.removesuffix(suffix)
            break
    else:
        file_base = Path(file_name).stem  # fallback for other extensions

    if args.humann_depleted_fasta: # file provided has already been depleted by humann (eg in bb4 workflows)
        file_base = file_base.replace("_bowtie2_unaligned", "")

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

    log_file = os.path.join(output_dir, file_base + "_baqlava.log")
    logging.basicConfig(
     filename=log_file,
     level=logging.INFO,
     format= '[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
     datefmt='%H:%M:%S',
     force=True)

    # Record versions of BAQLaVa, its sub-tools, and the databases in use.
    log_versions(humann_version)


    ########### File processing module ###########

    if args.bypass_bacterial_depletion: # bypassing bacterial depletion
        # NO UPSTREAM BACTERIAL DEPLETION

            workflow.add_task(
                "scp [depends[0]] [targets[0]]",
                depends = [args.input],
                targets = [str(tempdir / f"{file_base}_processed.fa")],
                name = "Upstream file processing")

    if args.humann_depleted_fasta: # file provided has already been depleted by humann (eg in bb4 workflows)
        workflow.add_task(
            "python [len_adj] [depends[0]] [targets[0]]",
            depends = [args.input],
            targets = [str(tempdir / f"{file_base}_processed.fa")],
            len_adj = os.path.abspath(args.lengthadjust),
            name = "Formatting bacterially depleted FASTA file")

    elif not args.bypass_bacterial_depletion: # bacterial depletion
        # UPSTREAM BACTERIAL DEPLETION
        # BACTERIAL DEPLETION (first step). This step fails gracefully: if HUMAnN
        # does not complete, the wrapper falls back to the original (non-depleted)
        # input so that nucleotide & translated search still proceed. The wrapper
        # always produces the processed FASTA target.
        if not args.taxonomic_profile:
            # PROCEED WITHOUT METAPHLAN TAXONOMIC PROFILE
            depletion_depends = [args.input]
            tax_profile = "NA"
        else:
            # USE A METAPHLAN TAXONOMIC PROFILE TO AID BACTERIAL DEPLETION
            depletion_depends = [args.input, args.taxonomic_profile]
            tax_profile = os.path.abspath(args.taxonomic_profile)

        workflow.add_task(
            "python [bd_script] [depends[0]] [tempdir] [base] [threads] [len_adj] [targets[0]] [tax]",
            depends = depletion_depends,
            targets = [str(tempdir / f"{file_base}_processed.fa")],
            bd_script = os.path.abspath(args.bacterial_depletion_script),
            tempdir = str(tempdir),
            base = file_base,
            threads = args.threads,
            len_adj = os.path.abspath(args.lengthadjust),
            tax = tax_profile,
            name = "Running HUMAnN to deplete bacterial reads from file")

        workflow.add_task(
            "scp [depends[0]] [targets[0]]",
            depends = [str(tempdir / f"{file_base}_processed.fa")],
            targets = [str(output_dir / f"{file_base}_bacterial_depleted.fa")],
            name = "Saving bacterially depleted FASTA file")

    ########### BAQLaVa NUCLEOIDE & TRANSLATED SEARCH ###########

    if not args.bypass_nucleotide_search: # bypass nucleotide != True

        # FIRST RUN & CALCULATE AT 25% COVERAGE:
        workflow.add_task(
            "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --count-normalization 'Adjusted RPKs' --nucleotide-subject-coverage-threshold 25 [args[2]]",
            depends = [str(tempdir / f"{file_base}_processed.fa")],
            args = [baq_dir, str(f"{file_base}_nucleotide1"), args.humann_passthrough_parameters_nucleotide],
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
            "humann --input [depends[0]] --output [args[0]] --bypass-nucleotide-index --nucleotide-database [n_db] --id-mapping [idx] --threads [threads] --bypass-translated-search --output-basename [args[1]] --count-normalization 'Adjusted RPKs' --nucleotide-subject-coverage-threshold 50 --resume [args[2]]",
            depends = [str(tempdir / f"{file_base}_processed.fa"), str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv")],
            args = [baq_dir, str(f"{file_base}_nucleotide2"), args.humann_passthrough_parameters_nucleotide],
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
            "humann --input [depends[0]] --output [args[0]] --id-mapping [idx] --protein-database [p_db] --threads [threads] --bypass-nucleotide-search --output-basename [args[1]] --count-normalization 'Adjusted RPKs' --translated-subject-coverage-threshold 50 [args[2]]",
            depends = [str(tempdir / f"{file_base}_processed.fa")],
            args = [baq_dir, str(f"{file_base}_translated"), args.humann_passthrough_parameters_translated],
            targets = [str(baq_dir / f"{file_base}_translated_2_genefamilies.tsv")],
            p_db = os.path.abspath(args.protdb),
            idx = os.path.abspath(translated_idmapping),
            threads = args.threads,
            name = "Running BAQLaVa Translated Search")


    ########### RECONCILING GENEFAMILIES.TSV FILES BASED ON PROCESSES RUN ###########

    if not args.bypass_nucleotide_search and args.bypass_translated_search: # nucleotide only, no translated

        workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [args[0]] [depends[2]] [proteomelen] [targets[0]] [targets[1]] [args[0]]",
            script = args.reconcile_mapped_script,
            mode = "1",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv"), str(baq_dir / f"{file_base}_nucleotide_50_genefamilies.tsv"), args.input],
            targets = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt"), str(output_dir / f"{file_base}_tempfile_markers.txt")],
            name = "Making BAQLaVa Viral Profile")

    elif args.bypass_nucleotide_search and not args.bypass_translated_search: # translated only, no nucleotide

        workflow.add_task("python [script] [mode] [args[0]] [args[0]] [depends[0]] [args[0]] [proteomelen] [targets[0]] [args[0]] [targets[1]]",
            script = args.reconcile_mapped_script,
            mode = "2",
            proteomelen = args.proteome_length,
            args = ["NA"],
            depends = [str(baq_dir / f"{file_base}_translated_2_genefamilies.tsv")],
            targets = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt"), str(output_dir / f"{file_base}_tempfile_proteins.txt")],
            name = "Making BAQLaVa Viral Profile")

    elif not args.bypass_nucleotide_search and not args.bypass_translated_search: # nucleotide + translated (full workflow)

            workflow.add_task("python [script] [mode] [depends[0]] [depends[1]] [depends[2]] [depends[3]] [proteomelen] [targets[0]] [targets[1]] [targets[2]]",
            script = args.reconcile_mapped_script,
            mode = "3",
            proteomelen = args.proteome_length,
            depends = [str(baq_dir / f"{file_base}_nucleotide_25_genefamilies.tsv"), str(baq_dir / f"{file_base}_nucleotide_50_genefamilies.tsv"), str(baq_dir / f"{file_base}_translated_2_genefamilies.tsv"), args.input],
            targets = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt"), str(output_dir / f"{file_base}_tempfile_markers.txt"), str(output_dir / f"{file_base}_tempfile_proteins.txt")],
            name = "Making BAQLaVa Viral Profile")

    else:
        print("ERROR: Incorrect reconciliation task prompted")

    ########### FILE CLEANUP ###########

    if not args.keep_tempfiles: # delete extra baqlava files
        workflow.add_task(
        "rm -r [args[0]]",
        depends = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt")],
        args = [baq_dir],
        name = "File Cleanup")

    workflow.add_task(
    "rm -r [args[0]]",
    depends = [str(output_dir / f"{file_base}_BAQLaVa_profile.txt")],
    args = [tempdir],
    name = "File Cleanup")

    # Nucleotide and translated search are independent (they share no task
    # dependency), so run them in parallel. Count how many search branches will
    # actually run (0, 1, or 2) and set the number of local jobs accordingly,
    # while honoring a higher user-supplied --local-jobs value. Steps within a
    # branch, depletion, reconciliation, and cleanup remain ordered by their
    # task dependencies regardless of the job count.
    run_nucleotide = not args.bypass_nucleotide_search
    run_translated = not args.bypass_translated_search
    parallel_search_jobs = int(run_nucleotide) + int(run_translated)
    local_jobs = max(args.jobs, parallel_search_jobs, 1)

    if run_nucleotide and run_translated and local_jobs > 1:
        message = "Running nucleotide and translated search in parallel (local jobs: " + str(local_jobs) + ")."
    elif run_nucleotide:
        message = "Running nucleotide search only."
    elif run_translated:
        message = "Running translated search only."
    else:
        message = "No search selected (both nucleotide and translated search bypassed)."
    logger.info(message)
    print("\n" + message + "\n")

    workflow.go(jobs=local_jobs)


if __name__ == '__main__':
    main()
