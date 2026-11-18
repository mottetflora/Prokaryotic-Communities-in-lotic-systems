#!/bin/bash


# First, we create a file with ORFs annotated with a KO involved in Carbon cycle
# NR==FNR to read the first file, meaning KO_used.txt
# ref[$1]=1 to take the first item of each line and store it in a table named ref
# $2 considers the second column in f file
# in ref { print $1 } checks if an ORF in the f file is present in the ref object
# f files are outputs of KofamScan, tables with ORFs and associated KOs
# KO_used.txt is a file with the 392 KO targeted involved the Carbon cycle
# The output is a list of every ORFs that were annotated with the targetted KOs

for f in /Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/*.txt; do
    Sample=$(basename "$f" .txt)

    awk 'NR==FNR { ref[$1]=1; next } $2 in ref { print $1 }' \
        /Users/floramottet/Documents/Pro/Bio-info/Script_MetaG/KO_used.txt "$f" > "/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORFs_KO/${Sample}_ORFs.txt"
done


# For every ORF annotated with a targeted KO (in File)
# seqkit greps the nucleic sequence associated to every ORF identified in the prodigal output (in Sequence) and store it in the output Sample_KO_seq.txt
for File in /Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORFs_KO/*_ORFs.txt; do

    Sample=$(basename "$File" _ORFs.txt)

    Sequence=/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Prodigal/${Sample}/*_output_mRNA.fa

    seqkit grep -f "$File" "$Sequence" > "/Users/floramottet/Documents/Pro/Bio-info/RE_metag/KO_seq/${Sample}_KO_seq.fa"

done

# kraken2 assigns cleaned raw reads against a taxonomic database
# --db : used database
# --report : results in "report" format, summarizing, lighter than output results files
# --use-names : returns used names and not TaxIDs
# Help : https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

for file in /beegfs/data/fmottet/MetaG_RE/KO_seq/*.fa ; do
	Sample="${file##*/}"
	echo $Sample
	/beegfs/data/fmottet/anaconda3/bin/kraken2 --db /ceph/banques/kraken/k2_core_nt_20251015 \
        	--report /beegfs/data/fmottet/MetaG_RE/KO_Kraken/${Sample}_kraken2_report --output /beegfs/data/fmottet/MetaG_RE/KO_Kraken/${Sample}_kraken2_results --use-names $file
done