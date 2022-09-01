## BAQLaVa V0

This is a very prelimiary version of the Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa). To run, your reads should be in a single fastq or fasta file (cat paired reads together into one file as needed). 

First, load the most recent version of biobakery workflows: 
  ```
  hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate
  ```
We reccomend depleting the fastq or fasta of potential bacterial reads before running baqlava. This can be done by running HUMAnN with standard parameters and skipping trasnlated search: 
  ```
  humann --input <FILE> --output <LOCATION> --bypass-trasnlated-search
  ```
The file at <LOCATION>/<FILE>_humann_temp/<FILE>_bowtie2_unaligned.fa has been depleted of bacterial reads via bowtie2 search to the chocophlan database. Format this file for BAQLaVa by running:
  ```
  python /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/remove_lengths_humann_bacterial_depletion.py
  ```
or by using the command line: 
  ```
  python /n/holystore01/LABS/huttenhower_lab/Users/jjensen/baqlava/run/remove_lengths_humann_bacterial_depletion.py
  ```
