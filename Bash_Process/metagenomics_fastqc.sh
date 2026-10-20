#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1  # core count
#SBATCH --partition=normal
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=flora.mottet@univ-lyon1.fr

# FastQC est un controle qualite des reads bruts.
# FastQC est deja present sur le cluster
# En premier, on met l'input qui doit etre dezippe
# -o c'est l'output, il doit deja avoir ete cree et doit etre vide
# --java stipule le chemin pour acceder a java (sinon, il ne le trouve pas)

cd /beegfs/data/fmottet/MetaG_RE/20251119-AVITI3-ABJE15WBR-0D2A15E2-C6D6-11F0-9DCF-DB86D3CEC64B/

for file in *_R1.fastq.gz; do
	sample=${file#*dir/}
	sample=${sample%_R1.fastq.gz}
	echo ${sample}
	gunzip ${sample}_R1.fastq.gz
	gunzip ${sample}_R2.fastq.gz
	### Change folder for results
	/beegfs/data/fmottet/Soft/FastQC/fastqc *.fastq -o /beegfs/data/fmottet/MetaG_RE/FastQC/DB86D3CEC64B/ --java /beegfs/data/soft/jdk1.8.0_201/bin/java
	gzip ${sample}_R1.fastq
	gzip ${sample}_R2.fastq
done

cd /beegfs/data/fmottet/MetaG_RE/20251128-AVITI3-ABS95IMHW-32CC64FA-CCDD-11F0-8AAE-62027A91FA7D/

for file in *_R1.fastq.gz; do
        sample=${file#*dir/}
        sample=${sample%_R1.fastq.gz}
        echo ${sample}
        gunzip ${sample}_R1.fastq.gz
        gunzip ${sample}_R2.fastq.gz
        ### Change folder for results
        /beegfs/data/fmottet/Soft/FastQC/fastqc *.fastq -o /beegfs/data/fmottet/MetaG_RE/FastQC/62027A91FA7D --java /beegfs/data/soft/jdk1.8.0_201/bin/java
        gzip ${sample}_R1.fastq
        gzip ${sample}_R2.fastq
done

cd /beegfs/data/fmottet/MetaG_RE/20251128-AVITI3-ABS95IMHW-3B4D9EAA-CCDD-11F0-B2D3-AD0B7A91FA7D/

for file in *_R1.fastq.gz; do
        sample=${file#*dir/}
        sample=${sample%_R1.fastq.gz}
        echo ${sample}
        gunzip ${sample}_R1.fastq.gz
        gunzip ${sample}_R2.fastq.gz
        ### Change folder for results
        /beegfs/data/fmottet/Soft/FastQC/fastqc *.fastq -o /beegfs/data/fmottet/MetaG_RE/FastQC/AD0B7A91FA7D --java /beegfs/data/soft/jdk1.8.0_201/bin/java
        gzip ${sample}_R1.fastq
        gzip ${sample}_R2.fastq
done

cd /beegfs/data/fmottet/MetaG_RE/20251201-AVITI3-AC1F3HHQT-FC5E16CA-D01A-11F0-B7C5-D18ED46030E8/

for file in *_R1.fastq.gz; do
        sample=${file#*dir/}
        sample=${sample%_R1.fastq.gz}
        echo ${sample}
        gunzip ${sample}_R1.fastq.gz
        gunzip ${sample}_R2.fastq.gz
        ### Change folder for results
        /beegfs/data/fmottet/Soft/FastQC/fastqc *.fastq -o /beegfs/data/fmottet/MetaG_RE/FastQC/D18ED46030E8 --java /beegfs/data/soft/jdk1.8.0_201/bin/java
        gzip ${sample}_R1.fastq
        gzip ${sample}_R2.fastq
done

