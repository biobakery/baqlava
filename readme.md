## BAQLaVa V0

This is a very prelimiary version of the Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa). To run, your reads should be in a single fastq or fasta file (cat paired reads together into one file as needed). 

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
Next, run HUMAnN with the BAQLaVa databases:
  ```
  humann --input <FILE> --output <LOCATION> --bypass-nucleotide-index --nucleotide-database /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/nucleotide_database_smallGVD --id-mapping /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/additional_files/idmap3.txt --bypass-translated-search
  ```
Use HUMAnN options (e.g. --threads) as you would like!

Finally, run the BAQLaVa reconciliation script. Make sure to reference the exact location of the input file (e.g. use ./ if it is in the current directory):
  ```
  python /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/reconcile_mapped_reads_nucleotide_only.py <FILE_with_LOCATION>
  ```
!! IMPORTANT !! Right now we are only using nucleotide mapping, not translated search. If you would like to use translated search, talk to Jordan about the protein database and setting that can be applied with:
  ```
  --protein-database /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/protein_database/
  ```
