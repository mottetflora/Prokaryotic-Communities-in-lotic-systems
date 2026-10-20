#!/bin/bash

# Script optimised by Mylene Hugoni and her team
# Xmx8g : memory to use (here 8 GB)
# in1 : raw forward reads 
# in2 : raw reverse reads 
# out1 : cleaned forward reads
# out2 : cleaned reverse reads
# minlen : minimal length of kept reads
# qtrim : right and left trimming
# trimq : minimal PHRED score of kept reads 
# ktrim : trimming on the right 
# k : k-mers length for reference base
# mink : adapters of 11bp long tested on kmers references
# ref : adapters file 
# tpe : trimming of reads R1 and R2 at the same length
# tbo : trimming of adapters on paired reads
# Documentation : https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

export PATH=/beegfs/data/soft/jdk1.8.0_201/bin:$PATH

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S3_clean_R2.fastq.gz \
        minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S4_clean_R2.fastq.gz \
        minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-GAB-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-GAB-A-S5_clean_R2.fastq.gz \
        minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MAT-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MAT-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-TON-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-TON-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MOL-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MOL-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PON-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PON-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-MGA-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-MGA-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-PLA-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-PLA-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-BRE-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-BRE-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo


bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S2_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S2_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S2_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S2_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S3_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S3_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S3_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S3_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S4_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S4_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S4_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S4_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo

bbduk.sh -Xmx8g threads=40 in1=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S5_R1.fastq.gz in2=/beegfs/data/fmottet/MetaG_RE/Data/AAP38-SAU-A-S5_R2.fastq.gz \
        out1=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S5_clean_R1.fastq.gz out2=/beegfs/data/fmottet/MetaG_RE/Clean/AAP38-SAU-A-S5_clean_R2.fastq.gz \
	minlen=50 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref=/beegfs/data/fmottet/git/Data/DB/adapters.fa hdist=1 tpe tbo
