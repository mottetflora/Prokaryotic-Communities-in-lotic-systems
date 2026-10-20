#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1  # core count
#SBATCH --partition=normal
#SBATCH -t 200:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=flora.mottet@etu.univ-lyon1.fr

# Script optimise par Mylene et son equipe
# Xmx8g : memoire a utiliser (ici 8 GB)
# in1 : le fichier de reads bruts dans le sens forward
# in2 : le fichier de reads bruts dans le sens reverse
# out1 : le fichier de sortie clean dans le sens forward
# out2 : le fichier de sortie clean dans le sens reverse
# minlen : taille minimum de reads qu'on garde
# qtrim : trimming du cote droit et gauche (rl : right et left)
# trimq : le score PHRED minimum des reads qu'on accepte de garder (selon Mylene, il vaut mieux choisir de garder beaucoup donc avec un faible score et filtrer apres l'assemblage de contigs)
# ktrim : avec r on trimme a droite avec l on trimme a gauche
# k : taille des k-mers qu'on va repertorier pour s'en servir de reference pour tester les adaptateurs
# mink : on testera les adaptateurs sur les k-mers de taille mink a k
# ref : fichiers avec les adaptateurs qu'on veut tester (donne par Mylene)
# hdist :
# tpe : trimming des reads R1 et R2 a la meme taille
# tbo : trimming des adapteurs sur les reads paires overlappes
# L'aide : https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

export PATH=/beegfs/data/soft/jdk1.8.0_201/bin:$PATH

# Pour les donnees de metagenomique

# Pour Gabriel
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


# Pour Mattant
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


# Pour Tonkin
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


# Pour Molottes
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


# Pour Ponton
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


# Pour Morte de Gabriel
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


# Pour Plaine
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


# Pour Morte de Brégnier
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


# Pour Saugey
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
