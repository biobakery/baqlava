## BAQLaVa V0.5

Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa) is made to be run on the FASRC severs with existing modules of HUMAnN & ANADAMA.

## Download & Installation:

To install BAQLaVa V0.5 on Harvard FASRC servers pleaase be using python 3.10, then run: 

     hutlab load rocky8/humann3/3.9-devel
     hutlab load rocky8/anadama2/0.10.0-devel
     git clone https://github.com/biobakery/baqlava.git
     python3 <DIR>/setup.py install --user

To install BAQLaVa V0.5 outside of the Harvard FASRC ecosystem you can clone the BAQLaVa repository and install. You will also need the following dependancies:

    HUMAnN 3.0 (https://github.com/biobakery/humann)
    AnADAMA2 (https://github.com/biobakery/anadama2)

When running BAQLaVa v0.5, you have the following options:
     
     --bypass-bacterial-depletion: run your sample through viral profiling without removing bacterial reads first 
     --bypass-nucleotide-search: run only the translated search component of BAQLaVa
     --bypass-translated-search: run only the nucleotide search component of BAQLaVa
     --taxonomic-profile: If you sample has previously been profiled with MetaPhlAn, you can speed up the bacterial depletion step and provide a MetaPhlAn taxonomic profile to be used directly
     
## Database requirements
Usage: (BAQLaVa DB + utility files are automatically downloaded and installed during setup)
```
baqlava_database --download database baqlava-db .
```

BAQLaVa requires two input databases: 1) a bowtie2-formatted nucleotide sequence database and 2) aa DIAMOND-formatted protein sequence database. BAQLaVa also needs two large reference files. These can be downloaded from https://huttenhower.sph.harvard.edu/baqlava-db/:
 
    https://huttenhower.sph.harvard.edu/baqlava-db/BAQLaVa.v1.0.hosted.tar.gz
   
Download the package above and unpack it. You will have a directory containing three subdirectories, each which need to be placed in a specific location under the install PATH for BAQLaVa:
    
    wget -P <LOCATION_TO_DOWNLOAD> https://g-227ca.190ebd.75bc.data.globus.org/baqlava-db/BAQLaVa.V0.5.tar.gz
    tar -zxvf BAQLaVa.V0.5.tar.gz
    mv BAQLaVa.V0.5/data/ <PATH>/baqlava/baqlava/.
    mv BAQLaVa.V0.5/utility_files/idmap_protein.txt <PATH>/baqlava/baqlava/utility_files/.
    mv BAQLaVa.V0.5/utility_files/nucleotide_marker_reference.txt <PATH>/baqlava/baqlava/utility_files/.
    mv BAQLaVa.V0.5/utility_files/translated_protein_reference.txt <PATH>/baqlava/baqlava/utility_files/.
    mv BAQLaVa.V0.5/utility_files/BAQLaVa_genomes_reference.txt <PATH>/baqlava/baqlava/utility_files/.

If you chose other locations to store this data, update the config file located at baqlava/baqlava/configs/baqlava.cfg so that the new locations for the protein and nucleotide databases and reference files are reflected. Otherwise, you do not need to update the config file, as it will assume the files are in the locations specified here.

## Input Data

To run, your reads should be in a single fastq or fasta file (cat paired end reads together into one file as needed). 

We reccomend using the depletion step to remove the fastq or fasta of potential bacterial reads before running baqlava. If you would like to bypass this step, you can include the flag: 
  ```
  --bypass-bacterial-depletion
  ```

## Running BAQLaVa

Finally, run BAQLaVa on your data!
```
baqlava -i <FILE> -o <OUTPUT_DIRECTORY>
```
BAQLaVa produces the following data products:

     FILE_BAQLaVa_profile.txt
     FILE_bacterial_depleted.fa
     FILE_tempfile_markers.txt
     FILE_tempfile_proteins.txt

FILE_BAQLaVa_profile.txt is the BAQLaVa viral profile. The other three files produced are not reqired for any further viral analysis but are produced to aid in future research. If you do not plan to utilize them, they can be discarded. 
     
     FILE_bacterial_depleted.fa:   The intermediate fasta file that has been bacterially depleted and only contains viral reads and other microbial dark matter reads.
     FILE_tempfile_markers.txt:    A list of all nucleotide markers mapped to, at what abundances, and which VGBs the markers belong to.
     FILE_tempfile_proteins.txt:   A list of all proteins mapped to, at what abundances, and which VGBs the proteins belong to.

If you'd like to access the full length BAQLaVa genomes, they are available for download here:

    https://g-227ca.190ebd.75bc.data.globus.org/baqlava-db/BAQLaVa.V0.5.raw_databases.tar.gz
     

