#!/usr/bin/env python
"""
BAQLaVa: bacterial depletion wrapper.

Runs HUMAnN to deplete bacterial reads from the input file, then formats the
resulting unaligned reads into the processed FASTA consumed by the downstream
nucleotide and translated search.

Bacterial depletion is the first (non-optional) BAQLaVa step, but it must fail
gracefully: if HUMAnN does not complete successfully, BAQLaVa falls back to the
original (non-depleted) input so that nucleotide and translated search can still
proceed. This script therefore always produces the processed-FASTA target and
exits 0, unless the fallback itself cannot be written.

Usage:
    run_bacterial_depletion.py <input> <tempdir> <file_base> <threads> \
        <lengthadjust_script> <processed_target> [taxonomic_profile]
"""

import os
import sys
import gzip
import shutil
import subprocess


def run(cmd):
    """Run a command, returning True on success (exit 0), False otherwise."""
    print("\nRunning: " + " ".join(str(c) for c in cmd) + "\n")
    try:
        return subprocess.run(cmd).returncode == 0
    except (FileNotFoundError, OSError) as e:
        sys.stderr.write("WARNING: Unable to run command (" + str(e) + "): "
                         + " ".join(str(c) for c in cmd) + "\n")
        return False


def copy_input_as_fallback(input_file, target):
    """Copy the original input to the processed-FASTA target (decompressing if
    the input is gzipped) so downstream search can run on the undepleted reads."""
    if input_file.endswith(".gz"):
        with gzip.open(input_file, "rt") as src, open(target, "wt") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(input_file, target)


def main():
    input_file = sys.argv[1]
    tempdir = sys.argv[2]
    file_base = sys.argv[3]
    threads = sys.argv[4]
    lengthadjust = sys.argv[5]
    target = sys.argv[6]
    taxonomic_profile = None
    if len(sys.argv) > 7 and sys.argv[7] not in ("", "NA", "False"):
        taxonomic_profile = sys.argv[7]

    unaligned = os.path.join(
        tempdir, file_base + "_humann_temp", file_base + "_bowtie2_unaligned.fa")

    humann_cmd = ["humann", "--input", input_file, "--output", tempdir,
                  "--bypass-translated-search", "--threads", str(threads)]
    if taxonomic_profile:
        humann_cmd += ["--taxonomic-profile", taxonomic_profile]

    depletion_ok = run(humann_cmd)

    if depletion_ok and os.path.isfile(unaligned):
        # Format the HUMAnN unaligned reads into the processed FASTA. The
        # length-adjust script takes (input, output_target).
        if not run(["python", lengthadjust, unaligned, target]):
            depletion_ok = False
    else:
        depletion_ok = False

    if depletion_ok and os.path.isfile(target):
        print("\nBacterial depletion completed successfully.\n")
        return

    # Fallback: depletion failed -> proceed with the original (undepleted) input.
    sys.stderr.write(
        "\nWARNING: Bacterial depletion failed. Proceeding to nucleotide and "
        "translated search using the original (non-depleted) input file.\n\n")
    try:
        copy_input_as_fallback(input_file, target)
    except (EnvironmentError, OSError) as e:
        sys.exit("ERROR: Bacterial depletion failed and the fallback input "
                 "could not be written (" + str(e) + "). Cannot continue.")


if __name__ == "__main__":
    main()
