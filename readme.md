  

# ***ATTENTION***

  

  

Before opening a new issue here, please check the appropriate help channel on the [bioBakery Support Forum](https://forum.biobakery.org) and consider opening or commenting on a thread there.

  

  

## BAQLaVa User Manual v0.5

  

Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa) is made to be run on the FASRC severs with existing modules of HUMAnN & ANADAMA.

  

  

## Contents ##
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data](#Input-Data)
- [How to Run](#how-to-run)
- [Demo Run](#Demo)
- [Output Data](#Output-Data)
- [Contributions](#Contributions)
- [Option List](#option-list)


## Requirements:

- HUMAnN>=3.9 (https://github.com/biobakery/humann)
- AnADAMA2>=0.10.0 (https://github.com/biobakery/anadama2)
- Python>=3.10 (https://github.com/biobakery/anadama2)
  

## Installation:

#### Install Using PyPi
```
pip install baqlava
```  
#### Install Using Github
```
git clone https://github.com/biobakery/baqlava.git
cd baqlava/
pip install .
```


Please use `--user` flag if you do not have permission to install the tools as a root user.
```
pip install baqlava --user
```

  

## Database requirements

  

  

Usage: (BAQLaVa DB instructions will be shown at the end of the setup completion as a banner message)

  

```
baqlava_database --download database baqlava-db /path/to/install_database
```

The above command will download and install the required BAQLaVa DB for the run in the specified path.

BAQLaVa requires two input databases: (`Total size ~1.6GB`)

  

1) a bowtie2-formatted nucleotide sequence database and

  

2) a DIAMOND-formatted protein sequence database.

  

#### Demo Database available here:

Nucleotide Demo DB: https://github.com/biobakery/baqlava/tree/master/examples/BAQLaVa.V0.5.nucleotide

Protein Demo DB: https://github.com/biobakery/baqlava/tree/master/examples/BAQLaVa.V0.5.protein

  

**NOTE:** If you'd like to access the full length BAQLaVa genomes, they are available for download here:

  

https://g-227ca.190ebd.75bc.data.globus.org/baqlava-db/BAQLaVa.V0.5.raw_databases.tar.gz

  

  

## Input Data

  

To run, your reads should be in a single fastq or fasta file (cat paired end reads together into one file as needed).

  

We reccomend using the depletion step to remove the fastq or fasta of potential bacterial reads before running baqlava. If you would like to bypass this step, you can include the flag `--bypass-bacterial-depletion`

  

Demo Input File: https://github.com/biobakery/baqlava/blob/master/examples/BAQLaVa_demo.fq

  

## How to run:

  

#### Demo Run

  

Finally, run BAQLaVa on your data!

  

```
baqlava -i <FILE> -o <OUTPUT_DIRECTORY>
```

  

  

When running BAQLaVa v0.5, you have the following options:

  

```
--bypass-bacterial-depletion: run your sample through viral profiling without removing bacterial reads first
--bypass-nucleotide-search: run only the translated search component of BAQLaVa
--bypass-translated-search: run only the nucleotide search component of BAQLaVa
--taxonomic-profile: If you sample has previously been profiled with MetaPhlAn, you can speed up the bacterial depletion step and provide a MetaPhlAn taxonomic profile to be used directly
```

  

  

## Output Data

  

BAQLaVa produces the following data products:

  

```
FILE_BAQLaVa_profile.txt
FILE_bacterial_depleted.fa
FILE_tempfile_markers.txt
FILE_tempfile_proteins.txt
```

  

  

`FILE_BAQLaVa_profile.txt` is the BAQLaVa viral profile. The other three files produced are not reqired for any further viral analysis but are produced to aid in future research. If you do not plan to utilize them, they can be discarded.

  

  

`FILE_bacterial_depleted.fa`: The intermediate fasta file that has been bacterially depleted and only contains viral reads and other microbial dark matter reads.

  

  

`FILE_tempfile_markers.txt`: A list of all nucleotide markers mapped to, at what abundances, and which VGBs the markers belong to.

  

  

`FILE_tempfile_proteins.txt`: A list of all proteins mapped to, at what abundances, and which VGBs the proteins belong to.

  

  

## Contributions ##

  

Thanks go to these wonderful people:

  

  

## Options List ##

  

All options can be accessed with `$ baqlava --help`.

  

```
usage: baqlava [-h] [--version] [--nucdb NUCDB] [--nucindex NUCINDEX] [--protindex PROTINDEX] [--threads THREADS] [--protdb PROTDB] [--lengthadjust LENGTHADJUST]
[--reconcile-mapped-script RECONCILE_MAPPED_SCRIPT] [--bypass-bacterial-depletion] [--bypass-nucleotide-search] [--bypass-translated-search]
[--taxonomic-profile TAXONOMIC_PROFILE] [--proteome-length PROTEOME_LENGTH] [--keep-tempfiles] -o OUTPUT [-i INPUT] [--config CONFIG] [--scripts SCRIPTS] [--tmp TMP]
[--local-jobs JOBS] [--grid-jobs GRID_JOBS] [--grid-tasks GRID_TASKS] [--grid GRID] [--grid-partition GRID_PARTITION] [--grid-benchmark {on,off}]  
[--grid-options GRID_OPTIONS] [--grid-submit-sleep GRID_SUBMIT_SLEEP] [--grid-environment GRID_ENVIRONMENT] [--grid-image GRID_IMAGE] [--grid-scratch GRID_SCRATCH]
[--grid-time-max GRID_TIME_MAX] [--grid-mem-max GRID_MEM_MAX] [--dry-run] [--skip-nothing] [--quit-early] [--until-task UNTIL_TASK] [--exclude-task EXCLUDE_TASK]
[--target TARGET] [--exclude-target EXCLUDE_TARGET] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

```