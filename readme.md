  

# ***ATTENTION***  

Before opening a new issue here, please check the appropriate help channel on the [bioBakery Support Forum](https://forum.biobakery.org) and consider opening or commenting on a thread there. 

For more in depth guidance on use, check out the BAQLaVa tutorial [here](https://github.com/biobakery/biobakery/wiki/BAQLaVa)! 

## BAQLaVa User Manual v0.5

Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa)

## Contents ##
- [Requirements](#requirements)
- [Installation](#installation)
- [How to Run](#Running-BAQLaVa)
- [Options](#Options-for-BAQLaVa-runs)
- [Output Data](#BAQLaVa-Output)
- [Creating a Merged Output from Individual Profiles](#Merging-BAQLaVa-Profiles)
- [Processing Paired Metagenomes + Metatranscriptomes](#Paired-MGX---MTX-Data)
- [BAQLaVa Viral Genomes](#BAQLaVa-Full-Length-Genomes)
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

BAQLaVa requires two input databases: (`Total size ~1.9GB`)

1) a bowtie2-formatted nucleotide sequence database (`~1.6GB`) and
  

2) a DIAMOND-formatted protein sequence database (`~0.3GB`).

These databases are large so you may want to specify a storage location for them outside of the default install path:

  
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

The main output is ```<FILENAME>_BAQLaVa_profile.txt``` which contains the viral profile. We will examine this in depth below. The other three files produced are not required for any further viral analysis but are produced to aid in future research. If you do not plan to utilize them, they can be discarded to save space.

```<FILENAME>_bacterial_depeled.fq``` is a copy of the input file, with any reads that mapped to a bacterial database having been removed, cutting down the number of reads within the file drastically. If you want to just look at possible viral reads for subsequent analysis in the future, this is a nice (much smaller!) file to work with.

```<FILENAME>_tempfile_markers.txt``` is a file that contains all markers that were mapped to in the nucleotide search step of BAQLaVa, the Viral Genome Bin (VGB) to which they belong, and their observed abundance.

```<FILENAME>_tempfile_proteins.txt``` is a file that contains all ORFs that were mapped to in the translated search step of BAQLaVa, the Viral Genome Bin (VGB) to which they belong, and their observed abundance.

## Options for BAQLaVa runs
When running BAQLaVa v0.5, you have the following options to modify the standard run:


```
--bypass-bacterial-depletion
```
 
Using the ```--bypass-bacterial-depletion``` flag will skip the first step of BAQLaVa which removes any bacterial reads from the input sample before proceeding to viral profiling. The bacterial depletion step is aimed at reducing false positive predictions. However, samples with particularly low abundance, along with metatranscriptomic samples, may require this mode if no bacterial reads are detected in the file. (Samples running through bacterial depletion but finding no bacterial reads will automatically fail. Rescue by rerunning with ```--bypass-bacterial-depletion```). 

```
--taxonomic-profile /PATH/TO/MetaPhlAn_profile
```

If a sample has already been profiled with MetaPhlAn previously and the bacteria present in a sample are therfore known, this flag can be used to speed up the bacterial depletion step. It will allow BAQLaVa to use the known profile to deplete bacterial reads rather than remapping entirely. Use the flag along with the path to a MetaPhlAn taxonomic profile. 


```
--bypass-nucleotide-search
--bypass-translated-search
```
These flags can be used to bypass either of the individual search steps of BAQLaVa if desired. We do not reccomend this as standard practice but may be useful for targeted research questions. You can always use the standard output stratified to the desired subset of information to obtain the same information.    


## BAQLaVa Output

The standard BAQLaVa profile output (```<FILENAME>_BAQLaVa_profile.txt```) looks like this for the demo data run above:
```
BAQLaVa VGB           BAQLaVa_demo_Abundance-RPKs    Reference Species      Taxonomy	                                                                                                                                                                  Other ICTV Genomes in VGB
VGB_1593              43.5046709396                  VGB_1593               r__duplodnaviria;k__heunggongvirae;p__uroviricota;c__caudoviricetes;o__crassvirales;f__crassvirales_unclassified;g__crassvirales_unclassified	
VGB_1593|nucleotide   13.184010221779264             VGB_1593               r__duplodnaviria;k__heunggongvirae;p__uroviricota;c__caudoviricetes;o__crassvirales;f__crassvirales_unclassified;g__crassvirales_unclassified	
VGB_1593|translated   43.5046709396                  VGB_1593               r__duplodnaviria;k__heunggongvirae;p__uroviricota;c__caudoviricetes;o__crassvirales;f__crassvirales_unclassified;g__crassvirales_unclassified	
VGB_49585             17.38481941893226              VGB_49585              r__duplodnaviria;k__heunggongvirae;p__uroviricota;c__caudoviricetes;o__caudoviricetes_unclassified;f__caudoviricetes_unclassified;g__caudoviricetes_unclassified	
VGB_49585|nucleotide  13.127554139030849             VGB_49585              r__duplodnaviria;k__heunggongvirae;p__uroviricota;c__caudoviricetes;o__caudoviricetes_unclassified;f__caudoviricetes_unclassified;g__caudoviricetes_unclassified	
VGB_49585|translated  17.38481941893226              VGB_49585              r__duplodnaviria;k__heunggongvirae;p__uroviricota;c__caudoviricetes;o__caudoviricetes_unclassified;f__caudoviricetes_unclassified;g__caudoviricetes_unclassified	
VGB_6438              39.447859546400004             VGB_6438               r__viruses_unclassified;k__viruses_unclassified;p__viruses_unclassified;c__viruses_unclassified;o__viruses_unclassified;f__viruses_unclassified;g__viruses_unclassified	
VGB_6438|translated   39.447859546400004             VGB_6438               r__viruses_unclassified;k__viruses_unclassified;p__viruses_unclassified;c__viruses_unclassified;o__viruses_unclassified;f__viruses_unclassified;g__viruses_unclassified

```

The BAQLaVa column ```BAQLaVa VGB``` names which Viral Genome Bin (VGB) was detected. Each VGB will have 2 or 3 lines. The first line will be the total BAQLaVa viral abundance of the VGB. The lines following the total abundance will show the abundance detected from each mapping step (```|nucleotide```, ```|translated```, or both). The total abundance is the max of the two individual abundances observed. The second column shows the individual or total abundance.
When Segment_Group is shown rather than a VGB, this indicates the viral species has a segmented genome, and as such represents a group of VGBs observed, but a single species bin.

The third column gives a reference species when the VGB contains any species formally recognized by ICTV. When the VGB does not contain a formally recognized species, the VGB name is given as the reference species. Here, no VGBs were identified which have an ICTV name.

The fourth column gives the full taxonomic lineage for the VGB. All lineages are ICTV-resolved and based on the known or predicted taxonomy of all genomes in the VGB. 

The final column will only have text when VGBs containing ICTV species genomes are identified (see column three for more info). When a VGB contains any ICTV species, the names of all ICTV-recognized species which are found in that VGB will be given in the fifth column.

## BAQLaVa Nuleotide & Translated Search: When to use individual vs. total abundances
We've tuned each of the individual seach steps (nucleotide & translated) to each serve a purpose, with nucleotide search being more specific, and translated search being more sensitive in the face of the growing amount of phage genomes being found that as a field we are far from having comprehensively captured or characterized. So, if the research question you have is something very focused, you may want to just look at the nucleotide subset lines. For general profiling, we recommend using the lines with the total viral abundance. 



## Merging BAQLaVa Profiles

BAQLaVa contains an accessory utility which can be used when mutliple samples have been profiled at a time. All individual sample BAQLaVa profiles should be located within a single directory, and then they can be merged with:
```
baqlava_join_tables /PATH/TO/DIRECTORY/ <FILE_PREFIX>
```

This utility will make a merged BAQLaVa profile of the total viral abundances across all samples, which will be named ```<FILE_PREFIX>_BAQLaVa_VGB_table.tsv```.


## Paired MGX - MTX Data

Because BAQLaVa works with both metagenomic and metatranscriptomic data, you may want to profile paired MGX-MTX data. To do so, we suggest the following workflow: 
1. QC both samples, making sure to remove rRNA reads from the MTX sample
2. Run MetaPhlan to carry out bacterial profiling on the MGX sample
3. Provide BAQLaVa the produced bacterial taxonomic profile for both the MTX and MGX sample (```--taxonomic-profile``` (see section above)). Doing this will deplete the same set of bacteria from both MGX and MTX samples.  


## BAQLaVa Full Length Genomes

The databases BAQLaVa uses in profiling are comprised of only VGB-specific genetic regions (markers) and ORFs. If you'd like to use the full genomes for the viruses within VGBs for further analysis, you can download them here:

https://g-227ca.190ebd.75bc.data.globus.org/baqlava-db/BAQLaVa.V0.5.raw_databases.tar.gz


  

## Contributions ##

  

Thanks go to these wonderful people:

Jordan Jensen, Eric Franzosa, Sagun Maharjan, Philipp Munch, Lea Wang, Curtis Huttenhower

  

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
