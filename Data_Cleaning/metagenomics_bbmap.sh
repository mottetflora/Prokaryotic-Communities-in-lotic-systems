#!/bin/bash

# bbmap is here used to map reads on predicted ORFs to then calculate ORFs cover
# ref with mRNA sequences ORFs predicted previously using prodigal
# in1 : raw forward reads 
# in2 : raw reverse reads 
# out : reads mapped against ORFs
# covstats : file of number of reads mapped against each ORF
# threads : number of threads to use
# nodisk : t to write nothing to on disk but outputs
# Documentation : https://bbmap.org/tools/bbmap 

DIR_RESULT="/beegfs/data/fmottet/MetaG_RE/ORF_Cover"

for folder in /beegfs/data/fmottet/MetaG_RE/Prodigal/*; do
	echo $folder
	sample="${folder##*/}"
	/beegfs/data/fmottet/Soft/bbmap/bbmap.sh ref=${folder}/${sample}_output_mRNA.fa  \
			in=/beegfs/data/fmottet/MetaG_RE/Clean/${sample}_clean_R1.fastq.gz  \
			in2=/beegfs/data/fmottet/MetaG_RE/Clean/${sample}_clean_R2.fastq.gz  \
                	out=${DIR_RESULT}/${sample}_prodigal_mapped.sam.gz \
                	covstats=${DIR_RESULT}/${sample}_prodigal_covstats.txt \
                	threads=16 \
                	nodisk=t
done
