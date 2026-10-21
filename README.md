# Prokaryotic-Communities-in-lotic-systems
Metabarcoding and metagenomic scripts for benthic prokaryotic communities structure and carbon cycling functioning

## Taxonomic structure of benthic procaryotic communities: Metabarcoding rRNA 16S data cleaning
1) Run in a bash terminal metabarcoding16S_bioinformatic_process.sh
2) Using R, run metabarcoding16S_cleaning_routine.R

## Functional profiles structure of benthic procaryotic structure: Metagenomics data cleaning
Steps 1 to 4 and 6 to 8 can be run in a bash terminal
1) Run metagenomics_fastqc.sh to check reads quality
2) Run metagenomics_bbduk.sh to remove biomolecular adapters
3) Run metagenomics_megahit.sh to assemble reads in contigs
4) Run metagenomics_quast.sh to check contigs quality
5) In a python environment, run metagenomics_contigs_filtration.py to filter contigs longer than 1 000 bp
6) Run metagenomics_prodigal.sh to perform prokaryotic Open Reading Frames prediction
7) Run metagenomics_bbmap.sh to map reads against ORFs (needed for ORFs cover)
8) Run metagenomics_kofamscan.sh to assign functions (KEGG Orthologs) to ORFs

## Statistical analyses
Using Load_Data.R script, create a R object with all data needed for statistical analyses presented in this study

All statistical analyses were performed using Statistical_Analysis_RE_paper.R

Metadata used in this study can be found in "Metadata" folder

Raw data can be found here: