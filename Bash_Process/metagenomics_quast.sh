#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1  # core count
#SBATCH --partition=normal
#SBATCH -t 2:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=flora.mottet@etu.univ-lyon1.fr

source /beegfs/data/fmottet/mypython/bin/activate

# Quast est un contrôle qualité des reads assemblés en contigs
# en premier on met l'input
# -o pour l'output

## Quast Megahit Co-assembly : 
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_all/final.contigs.fa \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/All_meg_quast

## Gabriel
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-GAB-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-GAB-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-GAB-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-GAB-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-GAB-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-GAB-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-GAB-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-GAB-A-S5

## Mattant
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MAT-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MAT-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MAT-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MAT-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MAT-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MAT-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MAT-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MAT-A-S5

## Tonkin
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-TON-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-TON-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-TON-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-TON-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-TON-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-TON-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-TON-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-TON-A-S5

## Molottes
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MOL-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MOL-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MOL-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MOL-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MOL-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MOL-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MOL-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MOL-A-S5

## Ponton
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PON-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PON-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PON-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PON-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PON-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PON-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PON-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PON-A-S5

## Morte de Gabriel
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MGA-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MGA-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MGA-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MGA-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MGA-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MGA-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-MGA-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-MGA-A-S5

## Plaine
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PLA-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PLA-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PLA-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PLA-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PLA-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PLA-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-PLA-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-PLA-A-S5

## Morte de Brégnier
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-BRE-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-BRE-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-BRE-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-BRE-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-BRE-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-BRE-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-BRE-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-BRE-A-S5

## Saugey
/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-SAU-A-S2/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-SAU-A-S2

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-SAU-A-S3/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-SAU-A-S3

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-SAU-A-S4/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-SAU-A-S4

/beegfs/data/fmottet/Soft/quast-5.3.0/quast.py /beegfs/data/fmottet/MetaG_RE/Megahit_results/AAP38-SAU-A-S5/contigs.fasta \
        -o /beegfs/data/fmottet/MetaG_RE/Quast/Megahit/AAP38-SAU-A-S5
