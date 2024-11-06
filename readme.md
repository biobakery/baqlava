## BAQLaVa V0.4

Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa) is made to be run on the FASRC severs with existing modules of HUMAnN & ANADAMA.

## Download & Installation:

To download BAQLaVa V0.5, first please request for the repository to be made public briefly. You can then run: 

     hutlab load rocky8/humann3/3.9-devel
     hutlab load rocky8/anadama2/0.10.0-devel
     git clone https://github.com/biobakery/baqlava.git
     python3 <DIR>/setup.py install --user

BAQLaVa v0.5 options:
     
     --bypass-bacterial-depletion
     --bypass-nucleotide-search
     --bypass-translated-search
     --taxonomic-profile (provide a MetaPhlAn taxonomic profile to aid bacterial depletion)
     
## Database requirements

BAQLaVa requires two input databases: 1) a bowtie2-formatted nucleotide sequence database and 2) aa DIAMOND-formatted protein sequence database. BAQLaVa also needs two large reference files. These can be downloaded from https://huttenhower.sph.harvard.edu/baqlava-db/:
 
    https://huttenhower.sph.harvard.edu/baqlava-db/BAQLaVa.v1.0.hosted.tar.gz
   
Download the package above and unpack it. You will have a directory containing three subdirectories, each which need to be placed in a specific location:
    
    wget -P <LOCATION_TO_DOWNLOAD> https://g-227ca.190ebd.75bc.data.globus.org/baqlava-db/BAQLaVa.V0.5.tar.gz
    tar -zxvf BAQLaVa.V0.5.tar.gz
    mv BAQLaVa.V0.5/data/ <PATH>/baqlava/baqlava/.
    mv BAQLaVa.V0.5/utility_files/idmap_protein.txt <PATH>/baqlava/baqlava/utility_files/.
    mv BAQLaVa.V0.5/utility_files/nucleotide_marker_reference.txt <PATH>/baqlava/baqlava/utility_files/.
    mv BAQLaVa.V0.5/utility_files/translated_protein_reference.txt <PATH>/baqlava/baqlava/utility_files/.

If you chose other locations to store this data, update the config file located at baqlava/baqlava/configs/baqlava.cfg so that the new locations for the protein and nucleotide databases and reference files are reflected.

## Input Data

To run, your reads should be in a single fastq or fasta file (cat paired reads together into one file as needed). 
First, load the most recent version of biobakery workflows: 
  ```
  hutlab load rocky8/humann3/3.9-devel
  hutlab load rocky8/anadama2/0.10.0-devel
  ```
We reccomend depleting the fastq or fasta of potential bacterial reads before running baqlava. If you would like to bypass this step, you can include the flag: 
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
     
     FILE_bacterial_depleted.fa:   The intermediate fasta file that has been bacterially depleted and only contains viral and microbial dark matter reads.
     FILE_tempfile_markers.txt:    A list of all nucleotide markers mapped to, at what abundances, and which VGBs the markers belong to.
     FILE_tempfile_proteins.txt:   A list of all proteins mapped to, at what abundances, and which VGBs the proteins belong to.
     

