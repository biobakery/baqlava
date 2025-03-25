  

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
  
Please ensure that HUMAnN is working properly before initiating a BAQLaVa run. 

## Installation:


### Step 1, Option 1: Install Using Github
This option provides demo data directly downloaded to test your install with.
```
git clone https://github.com/biobakery/baqlava.git
cd baqlava/
pip install .
```
#### Test the BAQLaVa install:

Before continuing, we will test the install of BAQLaVa is functional with a small demo database and demo data provided with the install:
```
baqlava -i baqlava/examples/baqlava_demo.fq -o <PATH/TO/OUTPUT> --nucdb examples/BAQLaVa.V0.5.nucleotide/ --protdb examples/BAQLaVa.V0.5.protein/
```
In the test run, we specify the demo databases to be supplied rather than full BAQLaVa databases. 

### Step 1, Option 2: Install Using PyPi
This option requires the demo data to be downloaded directly if you want to test your installation before moving on to download the full BAQLaVa databases. 
```
pip install baqlava
```  

Please use `--user` flag if you do not have permission to install the tools as a root user.
```
pip install baqlava --user
```
#### Download Demo Databases & Test the BAQLaVa Install:

Download the small demo data & databases:

Demo Input File: https://github.com/biobakery/baqlava/blob/master/examples/BAQLaVa_demo.fq

Nucleotide Demo DB: https://github.com/biobakery/baqlava/tree/master/examples/BAQLaVa.V0.5.nucleotide

Protein Demo DB: https://github.com/biobakery/baqlava/tree/master/examples/BAQLaVa.V0.5.protein

Test the BAQLaVa install:
```
baqlava -i <PATH/TO/FILE>baqlava_demo.fq -o <PATH/TO/OUTPUT> --nucdb /PATH/TO/BAQLaVa.V0.5.nucleotide/ --protdb /PATH/TO/BAQLaVa.V0.5.protein/
```
In the test run, we specify the demo databases to be supplied rather than full BAQLaVa databases.



### Step 2: Download Full BAQLaVa Databases:

Installing BAQLaVa from above will deliver a functional version of BAQLaVa with access to small demo databases only. In order to run BAQLaVa, download and install the full databases. The instructions to do so are below, and will also be shown at the end of the install as a banner message upon setup completion. 

BAQLaVa requires two input databases: (`Total size ~1.6GB`)

1) a bowtie2-formatted nucleotide sequence database and
  

2) a DIAMOND-formatted protein sequence database.

These databases are large so you may want to specify a storage location for them outside of the default install path!

  
Download and install the required BAQLaVa databases with:
```
baqlava_database --download database baqlava-db /path/to/install_database
```
With this completed, your BAQLaVa install should be complete and you should be ready to process samples!



## Running BAQLaVa

  

To run, your reads should be in a single fastq or fasta file. If you have paired end reads, we reccomend using ```cat``` to join paired end reads together into one file. To run, you only need to supply BAQLaVa the input file and the name of an output location for it to save outputs to. It will create this directory in the location specified if it does not already exist.

```
baqlava -i <FILE> -o <OUTPUT_DIRECTORY>
```

Running BAQLaVa creates four output products:
```
<FILENAME>_BAQLaVa_profile.txt
<FILENAME>_bacterial_depeled.fq
<FILENAME>_tempfile_markers.txt
<FILENAME>_tempfile_proteins.txt
```

The main output is ```<FILENAME>_BAQLaVa_profile.txt``` which contains the viral profile. We will examine this in depth below.

The additional files are products you may find helpful in further or downstream analysis of your viromes: 

```<FILENAME>_bacterial_depeled.fq``` is a copy of the input file, with any reads that mapped to a bacterial database having been removed, cutting down the number of reads within the file drastically. 

```<FILENAME>_tempfile_markers.txt``` is a file that contains all markers that were mapped to in the nucleotide search step of BAQLaVa, the Viral Genome Bin (VGB) to which they belong, and their observed abundance.

```<FILENAME>_tempfile_proteins.txt``` is a file that contains all ORFs that were mapped to in the translated search step of BAQLaVa, the Viral Genome Bin (VGB) to which they belong, and their observed abundance.

### Options for BAQLaVa runs
When running BAQLaVa v0.5, you have the following options to modify the standard run:


```
--bypass-bacterial-depletion
```
 
Using the ```--bypass-bacterial-depletion``` flag will skip the first step of BAQLaVa to remove any bacterial reads from the input sample before proceeding to viral profiling. The bacterial depletion step is aimed at reducing false positive predictions. However, samples with particularly low abundance, along with metatranscriptomic samples, may require this mode if no bacterial reads are detected in the file. (Samples running through bacterial depletion but finding no bacterial reads will automatically fail. Rescue by rerunning with ```--bypass-bacterial-depletion```). 

```
--taxonomic-profile /PATH/TO/MetaPhlAn_profile
```

If a sample has already been profiled with MetaPhlAn previously and as such the bacteria present in a sample are known, this flag can be used to speed up the bacterial depletion step. It will allow BAQLaVa to use a known profile to deplete bacterial reads rather than remapping entirely. Use the flag along with the path to a MetaPhlAn taxonomic profile. 


```
--bypass-nucleotide-search
--bypass-translated-search
```
These flags can be used to bypass either of the individual search steps of BAQLaVa if desired. We do not reccomend this as standard practice but may be useful for targeted research questions.  




  

**NOTE:** If you'd like to access the full length BAQLaVa genomes, they are available for download here:

  

https://g-227ca.190ebd.75bc.data.globus.org/baqlava-db/BAQLaVa.V0.5.raw_databases.tar.gz

  

  

## Input Data

  

To run, your reads should be in a single fastq or fasta file (cat paired end reads together into one file as needed).

  

We reccomend using the depletion step to remove the fastq or fasta of potential bacterial reads before running baqlava. If you would like to bypass this step, you can include the flag `--bypass-bacterial-depletion`

  

  

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
