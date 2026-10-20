#!/bin/python3.11
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1  # core count
#SBATCH --partition=normal
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=flora.mottet@etu.univ-lyon1.fr

import sys
sys.path.append("/beegfs/data/fmottet/mypython/lib/python3.11/site-packages")

# Aide : https://biopython.org/docs/1.75/api/Bio.SeqIO.html

import Bio.SeqIO
def filtration(input, size, output) :
        records = Bio.SeqIO.parse(input, "fasta")
        contigs_filtres = list()
        for record in records :
                if len(record.seq) >= size :
                        contigs_filtres.append(record)
        Bio.SeqIO.write(contigs_filtres, output, "fasta")
        print(output)

## GAB
print("Megahit_AAP38-GAB-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-GAB-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-GAB-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-GAB-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-GAB-A-S5/filtered2500_contigs.fasta")

## MAT
print("Megahit_AAP38-MAT-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-MAT-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-MAT-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-MAT-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MAT-A-S5/filtered2500_contigs.fasta")


## TON
print("Megahit_AAP38-TON-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-TON-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-TON-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-TON-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-TON-A-S5/filtered2500_contigs.fasta")


## MOL
print("Megahit_AAP38-MOL-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-MOL-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-MOL-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-MOL-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MOL-A-S5/filtered2500_contigs.fasta")


## PON
print("Megahit_AAP38-PON-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-PON-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-PON-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-PON-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PON-A-S5/filtered2500_contigs.fasta")


## MGA
print("Megahit_AAP38-MGA-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-MGA-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-MGA-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-MGA-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-MGA-A-S5/filtered2500_contigs.fasta")


## PLA
print("Megahit_AAP38-PLA-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-PLA-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-PLA-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-PLA-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-PLA-A-S5/filtered2500_contigs.fasta")


## BRE
print("Megahit_AAP38-BRE-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-BRE-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-BRE-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-BRE-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-BRE-A-S5/filtered2500_contigs.fasta")


## SAU
print("Megahit_AAP38-SAU-A-S2")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S2/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S2/filtered2500_contigs.fasta")

print("Megahit_AAP38-SAU-A-S3")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S3/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S3/filtered2500_contigs.fasta")

print("Megahit_AAP38-SAU-A-S4")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S4/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S4/filtered2500_contigs.fasta")

print("Megahit_AAP38-SAU-A-S5")
filtration("/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S5/final.contigs.fa", 2500, "/beegfs/data/fmottet/MetaG_RE/Megahit_results/Megahit_AAP38-SAU-A-S5/filtered2500_contigs.fasta")
