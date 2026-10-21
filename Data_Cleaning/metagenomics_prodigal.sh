#!/bin/bash

# Prodigal is a proteins predictor
# -p : anonymous for metagenomes and genomes
# -i : filtered contigs as input
# -o : general output
# -a : output file with proteins sequences
# -d : output file with mRNA sequences
# Documentation : https://www.biocode.ltd/catalog-tooldata/prodigal && https://github.com/hyattpd/prodigal/wiki

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S2/AAP38-GAB-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S2/AAP38-GAB-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S2/AAP38-GAB-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S3/AAP38-GAB-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S3/AAP38-GAB-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S3/AAP38-GAB-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S4/AAP38-GAB-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S4/AAP38-GAB-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S4/AAP38-GAB-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S5/AAP38-GAB-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S5/AAP38-GAB-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-GAB-A-S5/AAP38-GAB-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S2/AAP38-MAT-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S2/AAP38-MAT-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S2/AAP38-MAT-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S3/AAP38-MAT-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S3/AAP38-MAT-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S3/AAP38-MAT-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S4/AAP38-MAT-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S4/AAP38-MAT-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S4/AAP38-MAT-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S5/AAP38-MAT-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S5/AAP38-MAT-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MAT-A-S5/AAP38-MAT-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S2/AAP38-TON-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S2/AAP38-TON-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S2/AAP38-TON-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S3/AAP38-TON-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S3/AAP38-TON-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S3/AAP38-TON-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S4/AAP38-TON-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S4/AAP38-TON-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S4/AAP38-TON-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S5/AAP38-TON-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S5/AAP38-TON-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-TON-A-S5/AAP38-TON-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S2/AAP38-MOL-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S2/AAP38-MOL-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S2/AAP38-MOL-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S3/AAP38-MOL-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S3/AAP38-MOL-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S3/AAP38-MOL-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S4/AAP38-MOL-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S4/AAP38-MOL-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S4/AAP38-MOL-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S5/AAP38-MOL-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S5/AAP38-MOL-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MOL-A-S5/AAP38-MOL-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S2/AAP38-PON-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S2/AAP38-PON-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S2/AAP38-PON-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S3/AAP38-PON-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S3/AAP38-PON-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S3/AAP38-PON-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S4/AAP38-PON-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S4/AAP38-PON-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S4/AAP38-PON-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S5/AAP38-PON-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S5/AAP38-PON-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PON-A-S5/AAP38-PON-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S2/AAP38-MGA-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S2/AAP38-MGA-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S2/AAP38-MGA-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S3/AAP38-MGA-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S3/AAP38-MGA-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S3/AAP38-MGA-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S4/AAP38-MGA-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S4/AAP38-MGA-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S4/AAP38-MGA-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S5/AAP38-MGA-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S5/AAP38-MGA-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-MGA-A-S5/AAP38-MGA-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S2/AAP38-PLA-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S2/AAP38-PLA-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S2/AAP38-PLA-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S3/AAP38-PLA-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S3/AAP38-PLA-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S3/AAP38-PLA-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S4/AAP38-PLA-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S4/AAP38-PLA-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S4/AAP38-PLA-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S5/AAP38-PLA-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S5/AAP38-PLA-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-PLA-A-S5/AAP38-PLA-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S2/AAP38-BRE-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S2/AAP38-BRE-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S2/AAP38-BRE-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S3/AAP38-BRE-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S3/AAP38-BRE-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S3/AAP38-BRE-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S4/AAP38-BRE-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S4/AAP38-BRE-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S4/AAP38-BRE-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S5/AAP38-BRE-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S5/AAP38-BRE-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-BRE-A-S5/AAP38-BRE-A-S5_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S2/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S2/AAP38-SAU-A-S2_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S2/AAP38-SAU-A-S2_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S2/AAP38-SAU-A-S2_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S3/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S3/AAP38-SAU-A-S3_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S3/AAP38-SAU-A-S3_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S3/AAP38-SAU-A-S3_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S4/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S4/AAP38-SAU-A-S4_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S4/AAP38-SAU-A-S4_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S4/AAP38-SAU-A-S4_output_mRNA.fa

prodigal -p meta -i /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S5/filtered1000_contigs.fasta \
    -o /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S5/AAP38-SAU-A-S5_output_file.gbk \
    -a /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S5/AAP38-SAU-A-S5_output_protein.faa \
    -d /beegfs/data/fmottet/MetaG_RE/Prodigal/AAP38-SAU-A-S5/AAP38-SAU-A-S5_output_mRNA.fa