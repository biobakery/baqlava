## BAQLaVa V0.3

This is a pre-publication version of the Bioinformatic Application of Quantification and Labeling of Viral Taxonomy (BAQLaVa), made to be run on the FASRC severs with existing modules of HUMAnN & ANADAMA. The final version of BAQLaVa may contain algorithmic differences.

## Download & Installation:

To download BAQLaVa V0.3, first please request for the repository to be made public briefly. You can then run: 

     hutlab load rocky8/humann3/3.9-devel
     hutlab load rocky8/anadama2/0.10.0-devel
     git clone https://github.com/biobakery/baqlava.git
     python3 <DIR>/setup.py install --user

To test the installation, run: 

     baqlava -i examples/baqlava_demo.fastq -o output/

BAQLaVa v0.3 options:
     
     --bypass-bacterial-depletion
     --bypass-nucleotide-search
     --bypass-translated-search
     
The BAQLaVa data product is located in a file called <SAMPLE_PREFIX>_BAQLaVa_profile.txt (so for the demo data, this is "baqlava_demo_BAQLaVa_profile.tsv"). The data will look like: 

    Virus                               baqlava_demo_Abundance-RPKs	Virus_metadata	          Database
    Punavirus|Punavirus P1              30.0470091397	               Escherichia phage P1	nucleotide
    Tequatrovirus|Tequatrovirus T4      15.0099591829	               Escherichia phage T4	nucleotide
    Tunavirus|unclassified              61.74688894626818		                              translated

## Database requirements

BAQLaVa requires two input databases: 1) a bowtie2-formatted nucleotide sequence database and 2) aa DIAMOND-formatted protein sequence database. BAQLaVa also needs two large reference files. These can be downloaded from https://huttenhower.sph.harvard.edu/baqlava-db/:

 
    https://huttenhower.sph.harvard.edu/baqlava-db/BAQLaVa.v0.3.tar.gz
   
Download the package above and unpack it. You will have a directory containing three subdirectories, each which need to be placed in a specific location:
    
    wget -P <LOCATION_TO_DOWNLOAD> https://huttenhower.sph.harvard.edu/baqlava-db/BAQLaVa.v0.3.tar.gz
    tar -zxvf BAQLaVa.v0.3.tar.gz
    mv baqlava_release_v0.3_hostedfiles/BAQLaVa.V0.2.nucleotide <PATH>/baqlava/baqlava/data/.
    mv baqlava_release_v0.3_hostedfiles/BAQLaVa.V0.2.protein <PATH>/baqlava/baqlava/data/.
    mv baqlava_release_v0.3_hostedfiles/utility_files/nucleotide_marker_reference.txt <PATH>/baqlava/baqlava/utility_files/.
    mv baqlava_release_v0.3_hostedfiles/utility_files/translated_protein_reference.txt <PATH>/baqlava/baqlava/utility_files/.

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



