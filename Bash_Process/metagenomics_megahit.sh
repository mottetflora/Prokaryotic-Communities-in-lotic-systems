#!/bin/bash

source /beegfs/data/fmottet/mypython/bin/activate # activated python environment needed

# Megahit assembles reads overlaping, creating contigs
# Documentation : https://www.metagenomics.wiki/tools/assembly/megahit
# -1 : forwards reads 
# -2 : reverse reads
# -o : output
# -t : threads number to use
# -m : proportion of available memory to use (here 50% of allowed memory) 
# --k-list : lengths of tested k-mers

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S5

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S2_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S2_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S2

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S3_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S3_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S3

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S4_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S4_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S4

/beegfs/data/fmottet/Soft/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S5_clean_R1.fastq.gz \
        -2 /beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S5_clean_R2.fastq.gz --k-list 33,55,77 -t 6 -m 0.5 -o /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S5
