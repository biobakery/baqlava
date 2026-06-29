
# History #
## v1.1.0 06-29-2026
* Added --genome-filtering options to remove contaminant genomes from BAQLaVa VGBs in reponse to users identifying a small number of poorly performing VGBs (suspected plasmid contamination).
* Added a genome_comparison_utility to allow visualization of differences between VGBs at different levels of filtering.

## v1.0.0 03-01-2026
* This version is the same as v0.5.0, but was updated to v1.0 status with the release of the BAQLaVa preprint.

## v0.5.0 01-21-2025
* Automatically installed baqlava databases before the setup
* Added baqlava_databases executable script to download baqlava DB
* Added --keep-tempfiles flag
* Adjusted low abundance coverage estimate to account for sample read length

## v0.4.9 05-31-2024

* Updated protein search to target VGB-specific proteins only
* Updated taxonomy naming (bug fix)
* Added marker & ORF coverage parameters
* Added marker & ORF percent of sites covered parameters
* Changed abundance calculation for segmented genomes

## v0.3.0 03-23-2024

* Final pre-publication of BAQLaVa
* Markerized nucleotide search
* UniClust90 protein families
* Taxonomy addded to all VGBs

## v0.2.0 09-18-2023 ##

* Updated to use of markerized database for nucleotide and protein search.

## v0.1.0 07-20-2023 ##

* Expanded nucleotide database & ICTV genomes only used in protein search.
