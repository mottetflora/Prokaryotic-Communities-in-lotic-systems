#!/bin/bash

######################################
#### Quality controle with FastQC ####
######################################
### Documentation : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
cd /beegfs/home/fmottet/Data_Metabar_RhonEco/

for file in *_R1.fastq.gz; do
	sample=${file#*dir/}
	sample=${sample%_R1.fastq.gz}
	echo ${sample}
	gunzip ${sample}_R1.fastq.gz
	gunzip ${sample}_R2.fastq.gz
	### Change folder for results
	/beegfs/data/fmottet/Soft/FastQC/fastqc *.fastq -o /beegfs/data/fmottet/All_Seasons_RhonEco/FastQC --java /beegfs/data/soft/jdk1.8.0_201/bin/java
	gzip ${sample}_R1.fastq
	gzip ${sample}_R2.fastq
done


######## Number of reads in raw data
cd /beegfs/home/fmottet/Data_Metabar_RhonEco/

for file in *R1.fastq.gz; do
        sample=${file#*dir/}
        echo $sample
        unzip_file=${sample%.gz}
        echo $unzip_file >> /beegfs/data/fmottet/Data_RhonEco/Nb_reads.txt
        cp /beegfs/home/fmottet/Data_Metabar_RhonEco/$sample /beegfs/data/fmottet/All_Seasons_RhonEco/Data_Metabar/${unzip_file}.gz
        gunzip /beegfs/data/fmottet/All_Seasons_RhonEco/Data_Metabar/${unzip_file}.gz
        /beegfs/data/fmottet/Soft/obitools/bin/obicount /beegfs/data/fmottet/All_Seasons_RhonEco/Data_Metabar/$unzip_file >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_Metabar/Nb_reads_raw_data.txt
        rm /beegfs/data/fmottet/All_Seasons_RhonEco/Data_Metabar/$unzip_file
done

###########################################
## Merge paired-end reads with Obitools4 ##
###########################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023 ; page 12
### obipairing command : --min-identity = minimum of identity after alignement between R1 and R2 --min-overlap = minimum of overlaping bases between R1 and R2 after alignement -F sample_R1.fastq -R sample_R2.fastq > sample_merged_reads.fastq

cd /beegfs/home/fmottet/Data_Metabar_RhonEco/

for file in *.fastq.gz; do
	sample=${file%_R*.fastq.gz} ### % means without _R*.fastq.gz
	sample1=${sample}_R1.fastq.gz
	sample2=${sample}_R2.fastq.gz
	gunzip ${sample1}
	to_merge_R1=${sample}_R1.fastq
	gunzip ${sample2}
	to_merge_R2=${sample}_R2.fastq
	/beegfs/data/fmottet/Soft/obitools/bin/obipairing --min-identity=0.8 --min-overlap=10 -F $to_merge_R1 \
		-R $to_merge_R2 > /beegfs/data/fmottet/All_Seasons_RhonEco/Merged_Reads/${sample}_merged.fastq
	echo /beegfs/data/fmottet/All_Seasons_RhonEco/Merged_Reads/${sample}_merged.fastq >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_merged.txt
	/beegfs/data/fmottet/Soft/obitools/bin/obicount ${sample}_merged.fastq >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_merged.txt
	gzip $to_merge_R1
	gzip $to_merge_R2
done

###########################################
############## Demultiplexing #############
###########################################
cd /beegfs/data/fmottet/All_Seasons_RhonEco/Merged_Reads/

for file in *-LTRP6_L001_merged.fastq; do
	abrev=${file%-LTRP6_L001_merged.fastq}
	echo $abrev
	sample=${file%_merged.fastq}
	echo $sample
	/beegfs/data/fmottet/Soft/obitools/bin/obimultiplex -e 3 -t /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Table/${abrev}_table.txt -u /beegfs/data/fmottet/All_Seasons_RhonEco/Demulti/Unidentified/${abrev}_unidentified.fastq $file > /beegfs/data/fmottet/All_Seasons_RhonEco/Demulti/${abrev}_demulti.fastq
	echo $file >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_demulti.txt
	/beegfs/data/fmottet/Soft/obitools/bin/obicount $file >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_demulti.txt
	echo /beegfs/data/fmottet/All_Seasons_RhonEco/Demulti/${abrev}_demulti.fastq >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_demulti.txt
	/beegfs/data/fmottet/Soft/obitools/bin/obicount /beegfs/data/fmottet/All_Seasons_RhonEco/Demulti/${abrev}_demulti.fastq >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_demulti.txt
done

#####################################
######### Concatenate files #########
#####################################
cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

cat /beegfs/data/fmottet/Soft/obitools/bin/obicount /beegfs/data/fmottet/All_Seasons_RhonEco/Demulti/*.fastq >> all_seasons_reads.fastq

### obicount is here used to control any loss during the concatenation step
echo "All reads" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_cat.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount --max-cpu=4 all_seasons_reads.fastq >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_cat.txt

#############################################
### Convert fastq in fasta with Obitools4 ###
#############################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### Following command with --fasta-output option converts our fastq file into a fasta file
/beegfs/data/fmottet/Soft/obitools/bin/obiconvert --fasta-output all_seasons_reads.fastq > all_seasons_reads.fasta

### obicount is here used to control any loss during the conversion step
echo "All Reads in all_reads.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount --max-cpu=4 all_seasons_reads.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

#############################################
## Annotate sequence length with Obitools4 ##
#############################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### obiannotate with --length option add sequence's length in the sequence header
/beegfs/data/fmottet/Soft/obitools/bin/obiannotate --length all_seasons_reads.fasta > all_seasons_reads_length.fasta

awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq = seq $1} END {if (seq) print length(seq)}' all_seasons_reads.fasta > /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/length_reads.txt
# with /^>/ we indicate that each sequence starts with >
# with if (seq) print length(seq) we print previous sequence's length when we meet a new sequence sequence
# with seq=""; next we reinitialize seq variable and go to the next line 
# with seq = seq $1 we add every nucleotides lines to seq variable to get complete sequence without having multiple lines
# with END {if (seq) print length(seq) we print the length of the last sequence in the file

##########################################
### Filter joined reads with Obitools4 ###
##########################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### obigrep is here used to remove all sequences "joined" (unaligned) after merging step
/beegfs/data/fmottet/Soft/obitools/bin/obigrep -p 'annotations.mode != "join"' all_seasons_reads_length.fasta > all_seasons_reads_aligned.fasta

echo "All Reads aligned in all_reads_aligned.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount all_seasons_reads_aligned.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

###########################################
### Filtration by length with Obitools4 ###
###########################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023 + Script from M2 lessons Script_M2_GE_OK_corrige.txt

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### obigrep is here used to remove all sequences shorter than -l and longer than -L. -s option is used to test the pattern '^[acgt]+$' against the sequence (to remove sequences with N in nucleotides)
/beegfs/data/fmottet/Soft/obitools/bin/obigrep -l 300 -L 600 -s '^[acgt]+$' all_seasons_reads_aligned.fasta > all_seasons_reads_filtred.fasta

echo "All Reads filtred in all_seasons_reads_filtred.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount all_seasons_reads_filtred.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

########################################
##### Dereplication with Obitools4 #####
########################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### obiuniq is here used to merge by -m sample all reads to dereplicate the reads and to remove singletons (--no-singleton)
/beegfs/data/fmottet/Soft/obitools/bin/obiuniq -m sample --no-singleton all_seasons_reads_filtred.fasta > all_seasons_reads_uniq.fasta

echo "All Reads dereplicated in all_seasons_reads_uniq.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount all_seasons_reads_uniq.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

########################################
## Shorten reads names with Obitools4 ##
########################################
### Documentation : OBITools-V4.pdf by Eric Coissac ; 17 january 2023

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### obiannotate is used to shorten sequence's header keeping only count and merged_sample informations with -k option
/beegfs/data/fmottet/Soft/obitools/bin/obiannotate -k count all_seasons_reads_uniq.fasta >  all_seasons_reads_uniq_short.fasta

echo "All Reads with shorter names in all_seasons_reads_uniq_short.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount all_seasons_reads_uniq_short.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

sed -e 's/1:N:0:[A-Z][A-Z]*//g' -e 's/{"count":/count=/g' -e 's/}/;/g' all_seasons_reads_uniq_short.fasta > all_seasons_reads_sumaclust_format.fasta

#####################################
##### Clustering with Sumaclust #####
#####################################

cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### Clustering by sumaclust
### sumaclust command : -p = use multiple threads -t = identity of 0.97 (meaning 97%) to cluster sequences -R = maximum ratio of a sequence to be considered as a variant of the most abundant one -e = sequence assigned to the cluster with the center sequence presenting the highest similarity score > threshold -F = output in fasta format output.fasta input.fasta

/beegfs/data/soft/bin/sumaclust -p 4 -t 0.97 -R 0.25 -e -F all_seasons_reads_sumaclust_cluster.fasta all_seasons_reads_sumaclust_format.fasta

/beegfs/data/fmottet/Soft/obitools/bin/obiannotate -k cluster -k cluster_center -k cluster_score -k cluster_weight all_seasons_reads_sumaclust_cluster.fasta > all_seasons_reads_sumaclust_rscript.fasta

echo "All Reads after sumaclust in all_seasons_reads_sumaclust_rscript.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount all_seasons_reads_sumaclust_rscript.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

####### Do a scp on local machine of all_seasons_reads_uniq.fasta and all_seasons_reads_sumaclust_rscript.fasta
####### Run Merge_Sumaclust_Info.R R script on local machine
####### Do a scp on cluster of all_seasons_reads_final_cluster.fasta 

echo "All Reads after sumaclust annotations in all_seasons_reads_final_cluster.fasta" >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt
/beegfs/data/fmottet/Soft/obitools/bin/obicount all_seasons_reads_final_cluster.fasta >> /beegfs/data/fmottet/All_Seasons_RhonEco/Data_RhonEco/Nb_reads_final_steps.txt

#####################################
####### Taxonomic Assignement #######
#####################################
cd /beegfs/data/fmottet/All_Seasons_RhonEco/Final_steps/

### Taxonomic assignment with VSearch
### vsearch command : --usearch_gloabl = input for global alignment search --db = database --lca_cutoff = 95% of matching hits for the last common ancestor --lcaout = output for last common ancestor results --id = sequence rejected if the identity is lower than 95% --uc = output for global results --uc_allhits = put all the matches higher than 95% of identity in the output file --strand = dereplicate "both" strands

/beegfs/data/fmottet/anaconda3/bin/vsearch --usearch_global all_seasons_reads_final_cluster.fasta --db /beegfs/data/fmottet/Data_RhonEco/DB/SILVA_138.2_SSURef_NR99_tax_silva_trunc_vsearch.fasta \
        --lca_cutoff 0.95 --lcaout all_seasons_reads_lca_95.txt \
        --id 0.95 --uc all_seasons_reads_assigned_95.uc --uc_allhits --strand both
