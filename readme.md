## BAQLaVa V0.2

This is a very prelimiary version (V0) of the Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa), made to be run on the FASRC severs with existing modules of HUMAnN & ANADAMA. 

## Download & Installation:

To download BAQLaVa V0.1, first please request for the repository to be made public briefly. You can then run: 

     hutlab load centos7/python3/anadama2/0.10.0-devel
     hutlab load centos7/python3/humann3/3.6-devel
     git clone https://github.com/biobakery/baqlava.git
     python3 <DIR>/setup.py install --user

To test the installation, run: 

     baqlava -i examples/baqlava_demo.fastq -o output/
     
The BAQLaVa data product is located in a file called <SAMPLE_PREFIX>_BAQLaVa_profile.txt (so for the demo data, this is "baqlava_demo_BAQLaVa_profile.tsv"). The data will look like: 

    Virus	                         baqlava_demo_Abundance-RPKs	Virus_metadata	          Database
    Punavirus|Punavirus P1	     30.0470091397	               Escherichia phage P1	nucleotide
    Tequatrovirus|Tequatrovirus T4 15.0099591829	               Escherichia phage T4	nucleotide
    Tunavirus|unclassified	     61.74688894626818		                              translated

## Database requirements

BAQLaVa requires two input databases: 1) a bowtie2-formatted nucleotide sequence database and 2) aa DIAMOND-formatted protein sequence database. The databases can be downloaded from https://huttenhower.sph.harvard.edu/baqlava-db/:

 
    https://huttenhower.sph.harvard.edu/baqlava-db/BAQLaVa.V0.1.nucleotide.tar.gz
    https://huttenhower.sph.harvard.edu/baqlava-db/BAQLaVa.V0.1.protein.tar.gz   
   
    


## Input Data

To run, your reads should be in a single fastq or fasta file (cat paired reads together into one file as needed). 

First, load the most recent version of biobakery workflows: 
  ```
  hutlab load centos7/python3/humann3/3.6-devel
  ```
We reccomend depleting the fastq or fasta of potential bacterial reads before running baqlava. This can be done by running HUMAnN with standard parameters and skipping trasnlated search: 
  ```
  humann --input <FILE> --output <LOCATION> --bypass-translated-search
  ```
The file at LOCATION/FILE_humann_temp/FILE_bowtie2_unaligned.fa has been depleted of bacterial reads via bowtie2 search to the chocophlan database. Format this file for BAQLaVa by running:
  ```
  python /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/remove_lengths_humann_bacterial_depletion.py
  ```
This script creates a new filewith the suffix .lengthremoved.fa. This can alternatively be done through the command line and a new name chosen for the file:
  ```
  sed -r 's/\|[0-9]+$//' <FILE> > <NEWFILE>
  ```
## Running BAQLaVa - Option 1: manual via command line:

BAQLaVa is run via HUMAnN with the BAQLaVa databases:
  ```
  humann --input <FILE> --output <LOCATION> --bypass-nucleotide-index --nucleotide-database /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/nucleotide_database_smallGVD --id-mapping /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/additional_files/idmap4.txt --protein-database /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/protein_database/
  ```
Use HUMAnN options (e.g. --threads) as you would like!

Finally, run the BAQLaVa reconciliation script. Make sure to reference the exact location of the input file (e.g. use ./ if it is in the current directory). This will create a new file with "baqlava_genefamilies.tsv" as the suffix. This contains nucleotide mapping at the Species or genome (for MAGs) level and protein mapping at the Genus level. 
  ```
  python /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/reconcile_mapped_reads_v0.2.py <FILE_with_LOCATION>
  ```
## Running BAQLaVa - Option 2: AnADAMa2 workflow:
  ```
  python run.py -i <FOLDER_WITH_ALL_FILES> -o <OUTPUT_DIRECTORY>
  ```
