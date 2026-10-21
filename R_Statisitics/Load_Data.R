##################################################################################################
########### Script to save metabarcoding and metagenomic data for statistical analysis ###########
##################################################################################################

################################
###### Requiered packages ######
################################
library(dplyr)
library(ape)
library(vegan)
library(dplyr)

#################################
###### Metabarcoding files ######
#################################
### Load environmental dataset
## Dataset
Metadata_file <- read.csv("/Users/floramottet/Documents/Pro/RhonEco/Bio-info/R_Studio_total/Fichier_metadata.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE, na.strings = c(c("", " "), "NA"))

## Creation of a column with only site information
Metadata_file$Site_saison <- Metadata_file$Site
Metadata_file$Site <- gsub("S[12345]-", "", Metadata_file$Site)
Metadata_file$Lone <- gsub("S[12345]-", "", Metadata_file$Lone)

## Remove samples not biologicaly consistent 
Metadata_with_S1 <- Metadata_file[!grepl("SER|MOR|ROS|VIE", Metadata_file$Lone),]
Metadata <- Metadata_with_S1[7:297,]

## Two files : The original one with quantity and quality data and the next file with only quantity data
Metadata_used <- Metadata[,c(13, 15:17, 19:22, 26:29, 31:32)]
rownames(Metadata_used) <- Metadata$Pcr


## Agregate rows for useful environmental data for PCoA
Metadata_used$Site_saison <- Metadata$Site_saison
Metadata_used$Typologie <- Metadata$Typologie
Metadata_used$Nom <- Metadata$Nom
Metadata_used$Saison <- Metadata$Saison
Metadata_used$Month_Typo <- Metadata$Month_typo
Metadata_used$Site <- Metadata$Site
Metadata_used$Lone <- Metadata$Lone

Metadata_used_agreg <- as.data.frame(Metadata_used %>%
  group_by(Site_saison) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),      
    across(where(is.character), ~ first(.x))           
  ))

rownames(Metadata_used_agreg) <- Metadata_used_agreg$Site_saison

## Load eDNA dataset

## Reads abundances
## Reads per Sample for each MOTUs
Fichier_abondance_sommees <- read.table("/Users/floramottet/Documents/Pro/Bio-info/Data/sumaclust_filtered_data_sum_NEW.txt", header = TRUE)
Fichier_abondance_sommees <- Fichier_abondance_sommees[,!grepl("SER|MOR|ROS|VIE", names(Fichier_abondance_sommees))]
sum(colSums(Fichier_abondance_sommees[, -c(1:23)])) # 6 031 388 reads kept

## Assignations table + Reads per MOTUs for each site
Table_ab_sum <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/Data/Communauté_finale_reads_sommes_NEW.csv", dec = ",")
Table_ab_sum[, 1] <- gsub("AAP38-", "", Table_ab_sum[, 1])
Table_ab_sum[, 1] <- gsub("_.*", "", Table_ab_sum[, 1])
rownames(Table_ab_sum) <- Table_ab_sum[, 1]
Table_ab_sum <- Table_ab_sum[, -1]
Table_ab_sum[] <- lapply(Table_ab_sum, as.numeric)
noms_colonnes_reads_ab_sum <- readLines("/Users/floramottet/Documents/Pro/Bio-info/Data/Noms_colonnes_ab_sum_NEW.txt")
colnames(Table_ab_sum) <- noms_colonnes_reads_ab_sum
rm(noms_colonnes_reads_ab_sum)
Table_ab_sum <- Table_ab_sum[!grepl("SER|MOR|ROS|VIE", rownames(Table_ab_sum)),]
Table_ab_sum <- Table_ab_sum[!grepl("S1", rownames(Table_ab_sum)),]

### Relative reads abundances
## Reads per Sample for each MOTUs
Fichier_abondance_relatives <- read.table("/Users/floramottet/Documents/Pro/Bio-info/Data/sumaclust_filtered_data_rra_NEW.txt", header = TRUE)
Fichier_abondance_relatives <- Fichier_abondance_relatives[,!grepl("SER|MOR|ROS|VIE", names(Fichier_abondance_relatives))]
rownames(Fichier_abondance_relatives) <- Fichier_abondance_relatives$id
sum(colSums(Fichier_abondance_relatives[, -c(1:23)]))

## Assignations table + Reads per MOTUs for each site
Table_ab_rel <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/Data/Communauté_finale_reads_rra_NEW.csv", dec = ",")
Table_ab_rel[, 1] <- gsub("AAP38-", "", Table_ab_rel[, 1])
Table_ab_rel[, 1] <- gsub("_.*", "", Table_ab_rel[, 1])
rownames(Table_ab_rel) <- Table_ab_rel[, 1]
Table_ab_rel <- Table_ab_rel[, -1]
Table_ab_rel[] <- lapply(Table_ab_rel, as.numeric)
noms_colonnes_reads_ab_rel <- readLines("/Users/floramottet/Documents/Pro/Bio-info/Data/Noms_colonnes_ab_rel_NEW.txt")
colnames(Table_ab_rel) <- noms_colonnes_reads_ab_rel
rm(noms_colonnes_reads_ab_rel)
Table_ab_rel <- Table_ab_rel[!grepl("SER|MOR|ROS|VIE", rownames(Table_ab_rel)),]
Table_ab_rel <- Table_ab_rel[!grepl("S1", rownames(Table_ab_rel)),]

### Taxonomic table
Assigned_OTU <- Fichier_abondance_relatives[,10:15]
Assigned_OTU[is.na(Assigned_OTU)] <- "Other"

### MOTUs table
MOTUs_ab_som <- Fichier_abondance_relatives[,23:313]
colnames(MOTUs_ab_som) <- Metadata$Pcr

##Remove Chloroplast and Mitochondria MOTUs

## Removing unassigned MOTUs at domain level
sum(Assigned_OTU$Domaine=="Other") # [1] 240 OTUs chloroplastiques
exclude_unassigned_domains <- c(which(Assigned_OTU$Domaine=="Other"))
Assigned_OTU <- Assigned_OTU[-exclude_unassigned_domains,]
MOTUs_ab_som <- MOTUs_ab_som[-exclude_unassigned_domains,]
Table_ab_sum <- Table_ab_sum[,-exclude_unassigned_domains]
rm(exclude_unassigned_domains)

## Removing Chloroplast                      
sum(Assigned_OTU$Ordre=="Chloroplast") # [1] 240 OTUs chloroplastiques
exclude_chloro <- c(which(Assigned_OTU$Ordre=="Chloroplast"))
Assigned_OTU <- Assigned_OTU[-exclude_chloro,]
MOTUs_ab_som <- MOTUs_ab_som[-exclude_chloro,]
Table_ab_sum <- Table_ab_sum[,-exclude_chloro]
rm(exclude_chloro)

## Removing Chloroplast
sum(Assigned_OTU$Famille == "Mitochondria") # [1] 12 OTUs mitochondriaux
exclude_mito <- c(which(Assigned_OTU$Famille == "Mitochondria"))
Assigned_OTU <- Assigned_OTU[-exclude_mito,]
MOTUs_ab_som <- MOTUs_ab_som[-exclude_mito,]
Table_ab_sum <- Table_ab_sum[,-exclude_mito]
rm(exclude_mito)

## Rarefaction 
Data_raref <- rrarefy(round(Table_ab_sum, 0), min(round(rowSums(Table_ab_sum), 0)))

Data_raref_per_site <- aggregate(Data_raref, by=list(substr(rownames(Data_raref), 4, 8)), mean)
rownames(Data_raref_per_site) <- Data_raref_per_site[,1]
Data_raref_per_site<-Data_raref_per_site[,-1]

Data_raref_per_site_per_season <- aggregate(Data_raref, by=list(substr(rownames(Data_raref), 1, 8)), mean)
rownames(Data_raref_per_site_per_season) <- Data_raref_per_site_per_season[,1]
Data_raref_per_site_per_season<-Data_raref_per_site_per_season[,-1]

Data_raref_per_lone_par_saison <- aggregate(Data_raref, by=list(substr(rownames(Data_raref), 1, 6)), mean)
rownames(Data_raref_per_lone_par_saison) <- Data_raref_per_lone_par_saison[,1]
Data_raref_per_lone_par_saison <- Data_raref_per_lone_par_saison[,-1]

Data_raref_per_lone <- aggregate(Data_raref, by=list(substr(rownames(Data_raref), 4,6)), mean)
rownames(Data_raref_per_lone) <- Data_raref_per_lone[,1]
Data_raref_per_lone <- Data_raref_per_lone[,-1]

## Remove unimportant files
rm(Fichier_abondance_relatives)
rm(Fichier_abondance_sommees)
rm(Metadata_with_S1)
rm(Metadata_file)


#################################
####### Metagenomic files #######
#################################
## Metadata
Metadata_metaG <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/RStudio/Metadata_MetaG_RE.csv", header = TRUE, sep = ";")
rownames(Metadata_metaG) <- Metadata_metaG$Samples

Pathway_KO_KEGG_list <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/Script_MetaG/KO_Pathway_Kegg.csv", sep = ";", header = T)

Metadata_DESeq <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/RStudio/Metadata_MetaG_RE_DESeq.csv", header = T, sep = ";")
rownames(Metadata_DESeq) <- Metadata_DESeq$Samples
Metadata_DESeq$Typo_court <- c(rep("Disconnected", 4), rep("Connected", 8), rep("Disconnected", 4), rep("Inter", 4),
                               rep("Disconnected", 4), rep("Inter", 4), rep("Disconnected", 4), rep("Connected", 4))

rownames(Metadata_DESeq) <- Metadata_DESeq$Samples
Metadata_DESeq <- Metadata_DESeq[sort(Metadata_DESeq$Samples),]
KO_Samples_Carbon_Cycle_raw_count <- KO_Samples_Carbon_Cycle_raw_count[,order(match(colnames(KO_Samples_Carbon_Cycle_raw_count), rownames(Metadata_DESeq)))]
KO_Samples_Carbon_Cycle_raw_count_sans_inter <- KO_Samples_Carbon_Cycle_raw_count[,-c(17:20, 25:28)]
rownames(KO_Samples_Carbon_Cycle_raw_count_sans_inter) <- KO_C_Cycle_list
Metadata_DESeq <- Metadata_DESeq[-c(17:20, 25:28),]
Metadata_DESeq$Typo_court <- c(rep("Disconnected", 4), rep("Connected", 8), rep("Disconnected", 4),
                               rep("Disconnected", 4), rep("Disconnected", 4), rep("Connected", 4))

Metadata_DESeq$Typo_court <- factor(Metadata_DESeq$Typo_court,
                                    levels = c("Disconnected", "Connected")
)


## Load Mapping Data
BRE_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-BRE-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(BRE_S2_mapping) <- sub(" .*", "", BRE_S2_mapping[[1]])
BRE_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-BRE-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(BRE_S3_mapping) <- sub(" .*", "", BRE_S3_mapping[[1]])
BRE_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-BRE-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(BRE_S4_mapping) <- sub(" .*", "", BRE_S4_mapping[[1]])
BRE_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-BRE-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(BRE_S5_mapping) <- sub(" .*", "", BRE_S5_mapping[[1]])

GAB_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-GAB-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(GAB_S2_mapping) <- sub(" .*", "", GAB_S2_mapping[[1]])
GAB_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-GAB-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(GAB_S3_mapping) <- sub(" .*", "", GAB_S3_mapping[[1]])
GAB_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-GAB-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(GAB_S4_mapping) <- sub(" .*", "", GAB_S4_mapping[[1]])
GAB_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-GAB-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(GAB_S5_mapping) <- sub(" .*", "", GAB_S5_mapping[[1]])

MAT_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MAT-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MAT_S2_mapping) <- sub(" .*", "", MAT_S2_mapping[[1]])
MAT_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MAT-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MAT_S3_mapping) <- sub(" .*", "", MAT_S3_mapping[[1]])
MAT_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MAT-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MAT_S4_mapping) <- sub(" .*", "", MAT_S4_mapping[[1]])
MAT_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MAT-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MAT_S5_mapping) <- sub(" .*", "", MAT_S5_mapping[[1]])

MGA_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MGA-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MGA_S2_mapping) <- sub(" .*", "", MGA_S2_mapping[[1]])
MGA_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MGA-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MGA_S3_mapping) <- sub(" .*", "", MGA_S3_mapping[[1]])
MGA_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MGA-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MGA_S4_mapping) <- sub(" .*", "", MGA_S4_mapping[[1]])
MGA_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MGA-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MGA_S5_mapping) <- sub(" .*", "", MGA_S5_mapping[[1]])

MOL_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MOL-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MOL_S2_mapping) <- sub(" .*", "", MOL_S2_mapping[[1]])
MOL_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MOL-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MOL_S3_mapping) <- sub(" .*", "", MOL_S3_mapping[[1]])
MOL_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MOL-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MOL_S4_mapping) <- sub(" .*", "", MOL_S4_mapping[[1]])
MOL_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-MOL-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(MOL_S5_mapping) <- sub(" .*", "", MOL_S5_mapping[[1]])

PLA_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PLA-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PLA_S2_mapping) <- sub(" .*", "", PLA_S2_mapping[[1]])
PLA_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PLA-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PLA_S3_mapping) <- sub(" .*", "", PLA_S3_mapping[[1]])
PLA_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PLA-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PLA_S4_mapping) <- sub(" .*", "", PLA_S4_mapping[[1]])
PLA_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PLA-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PLA_S5_mapping) <- sub(" .*", "", PLA_S5_mapping[[1]])

PON_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PON-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PON_S2_mapping) <- sub(" .*", "", PON_S2_mapping[[1]])
PON_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PON-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PON_S3_mapping) <- sub(" .*", "", PON_S3_mapping[[1]])
PON_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PON-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PON_S4_mapping) <- sub(" .*", "", PON_S4_mapping[[1]])
PON_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-PON-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(PON_S5_mapping) <- sub(" .*", "", PON_S5_mapping[[1]])

SAU_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-SAU-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(SAU_S2_mapping) <- sub(" .*", "", SAU_S2_mapping[[1]])
SAU_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-SAU-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(SAU_S3_mapping) <- sub(" .*", "", SAU_S3_mapping[[1]])
SAU_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-SAU-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(SAU_S4_mapping) <- sub(" .*", "", SAU_S4_mapping[[1]])
SAU_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-SAU-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(SAU_S5_mapping) <- sub(" .*", "", SAU_S5_mapping[[1]])

TON_S2_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-TON-A-S2_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(TON_S2_mapping) <- sub(" .*", "", TON_S2_mapping[[1]])
TON_S3_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-TON-A-S3_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(TON_S3_mapping) <- sub(" .*", "", TON_S3_mapping[[1]])
TON_S4_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-TON-A-S4_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(TON_S4_mapping) <- sub(" .*", "", TON_S4_mapping[[1]])
TON_S5_mapping <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/ORF_Cover/AAP38-TON-A-S5_prodigal_covstats.csv", sep = ";", header = TRUE)
rownames(TON_S5_mapping) <- sub(" .*", "", TON_S5_mapping[[1]])

## Load Assignation Data
BRE_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/BRE-A-S2_kofam.csv", sep = ";", header = FALSE)
BRE_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/BRE-A-S3_kofam.csv", sep = ";", header = FALSE)
BRE_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/BRE-A-S4_kofam.csv", sep = ";", header = FALSE)
BRE_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/BRE-A-S5_kofam.csv", sep = ";", header = FALSE)

GAB_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/GAB-A-S2_kofam.csv", sep = ";", header = FALSE)
GAB_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/GAB-A-S3_kofam.csv", sep = ";", header = FALSE)
GAB_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/GAB-A-S4_kofam.csv", sep = ";", header = FALSE)
GAB_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/GAB-A-S5_kofam.csv", sep = ";", header = FALSE)

MAT_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MAT-A-S2_kofam.csv", sep = ";", header = FALSE)
MAT_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MAT-A-S3_kofam.csv", sep = ";", header = FALSE)
MAT_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MAT-A-S4_kofam.csv", sep = ";", header = FALSE)
MAT_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MAT-A-S5_kofam.csv", sep = ";", header = FALSE)

MGA_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MGA-A-S2_kofam.csv", sep = ";", header = FALSE)
MGA_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MGA-A-S3_kofam.csv", sep = ";", header = FALSE)
MGA_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MGA-A-S4_kofam.csv", sep = ";", header = FALSE)
MGA_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MGA-A-S5_kofam.csv", sep = ";", header = FALSE)

MOL_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MOL-A-S2_kofam.csv", sep = ";", header = FALSE)
MOL_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MOL-A-S3_kofam.csv", sep = ";", header = FALSE)
MOL_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MOL-A-S4_kofam.csv", sep = ";", header = FALSE)
MOL_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/MOL-A-S5_kofam.csv", sep = ";", header = FALSE)

PLA_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PLA-A-S2_kofam.csv", sep = ";", header = FALSE)
PLA_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PLA-A-S3_kofam.csv", sep = ";", header = FALSE)
PLA_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PLA-A-S4_kofam.csv", sep = ";", header = FALSE)
PLA_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PLA-A-S5_kofam.csv", sep = ";", header = FALSE)

PON_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PON-A-S2_kofam.csv", sep = ";", header = FALSE)
PON_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PON-A-S3_kofam.csv", sep = ";", header = FALSE)
PON_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PON-A-S4_kofam.csv", sep = ";", header = FALSE)
PON_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/PON-A-S5_kofam.csv", sep = ";", header = FALSE)

SAU_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/SAU-A-S2_kofam.csv", sep = ";", header = FALSE)
SAU_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/SAU-A-S3_kofam.csv", sep = ";", header = FALSE)
SAU_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/SAU-A-S4_kofam.csv", sep = ";", header = FALSE)
SAU_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/SAU-A-S5_kofam.csv", sep = ";", header = FALSE)

TON_S2_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/TON-A-S2_kofam.csv", sep = ";", header = FALSE)
TON_S3_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/TON-A-S3_kofam.csv", sep = ";", header = FALSE)
TON_S4_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/TON-A-S4_kofam.csv", sep = ";", header = FALSE)
TON_S5_KO <- read.csv("/Users/floramottet/Documents/Pro/Bio-info/RE_metag/Kofam/TON-A-S5_kofam.csv", sep = ";", header = FALSE)

## Aggregate Mapping and KO datasets
BRE_S2_KO_merged <- aggregate(V2 ~ V1, data = BRE_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
BRE_S2_mapping$ORF <- rownames(BRE_S2_mapping)
BRE_S2_mapping_KO <- merge(BRE_S2_mapping, BRE_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(BRE_S2_mapping_KO) <- BRE_S2_mapping_KO$ORF
BRE_S2_mapping_KO[BRE_S2_mapping_KO == ""] <- NA
rm(BRE_S2_KO, BRE_S2_KO_merged, BRE_S2_mapping)

BRE_S3_KO_merged <- aggregate(V2 ~ V1, data = BRE_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
BRE_S3_mapping$ORF <- rownames(BRE_S3_mapping)
BRE_S3_mapping_KO <- merge(BRE_S3_mapping, BRE_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(BRE_S3_mapping_KO) <- BRE_S3_mapping_KO$ORF
BRE_S3_mapping_KO[BRE_S3_mapping_KO == ""] <- NA
rm(BRE_S3_KO, BRE_S3_KO_merged, BRE_S3_mapping)

BRE_S4_KO_merged <- aggregate(V2 ~ V1, data = BRE_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
BRE_S4_mapping$ORF <- rownames(BRE_S4_mapping)
BRE_S4_mapping_KO <- merge(BRE_S4_mapping, BRE_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(BRE_S4_mapping_KO) <- BRE_S4_mapping_KO$ORF
BRE_S4_mapping_KO[BRE_S4_mapping_KO == ""] <- NA
rm(BRE_S4_KO, BRE_S4_KO_merged, BRE_S4_mapping)

BRE_S5_KO_merged <- aggregate(V2 ~ V1, data = BRE_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
BRE_S5_mapping$ORF <- rownames(BRE_S5_mapping)
BRE_S5_mapping_KO <- merge(BRE_S5_mapping, BRE_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(BRE_S5_mapping_KO) <- BRE_S5_mapping_KO$ORF
BRE_S5_mapping_KO[BRE_S5_mapping_KO == ""] <- NA
rm(BRE_S5_KO, BRE_S5_KO_merged, BRE_S5_mapping)



GAB_S2_KO_merged <- aggregate(V2 ~ V1, data = GAB_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
GAB_S2_mapping$ORF <- rownames(GAB_S2_mapping)
GAB_S2_mapping_KO <- merge(GAB_S2_mapping, GAB_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(GAB_S2_mapping_KO) <- GAB_S2_mapping_KO$ORF
GAB_S2_mapping_KO[GAB_S2_mapping_KO == ""] <- NA
rm(GAB_S2_KO, GAB_S2_KO_merged, GAB_S2_mapping)

GAB_S3_KO_merged <- aggregate(V2 ~ V1, data = GAB_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
GAB_S3_mapping$ORF <- rownames(GAB_S3_mapping)
GAB_S3_mapping_KO <- merge(GAB_S3_mapping, GAB_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(GAB_S3_mapping_KO) <- GAB_S3_mapping_KO$ORF
GAB_S3_mapping_KO[GAB_S3_mapping_KO == ""] <- NA
rm(GAB_S3_KO, GAB_S3_KO_merged, GAB_S3_mapping)

GAB_S4_KO_merged <- aggregate(V2 ~ V1, data = GAB_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
GAB_S4_mapping$ORF <- rownames(GAB_S4_mapping)
GAB_S4_mapping_KO <- merge(GAB_S4_mapping, GAB_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(GAB_S4_mapping_KO) <- GAB_S4_mapping_KO$ORF
GAB_S4_mapping_KO[GAB_S4_mapping_KO == ""] <- NA
rm(GAB_S4_KO, GAB_S4_KO_merged, GAB_S4_mapping)

GAB_S5_KO_merged <- aggregate(V2 ~ V1, data = GAB_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
GAB_S5_mapping$ORF <- rownames(GAB_S5_mapping)
GAB_S5_mapping_KO <- merge(GAB_S5_mapping, GAB_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(GAB_S5_mapping_KO) <- GAB_S5_mapping_KO$ORF
GAB_S5_mapping_KO[GAB_S5_mapping_KO == ""] <- NA
rm(GAB_S5_KO, GAB_S5_KO_merged, GAB_S5_mapping)



MAT_S2_KO_merged <- aggregate(V2 ~ V1, data = MAT_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MAT_S2_mapping$ORF <- rownames(MAT_S2_mapping)
MAT_S2_mapping_KO <- merge(MAT_S2_mapping, MAT_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MAT_S2_mapping_KO) <- MAT_S2_mapping_KO$ORF
MAT_S2_mapping_KO[MAT_S2_mapping_KO == ""] <- NA
rm(MAT_S2_KO, MAT_S2_KO_merged, MAT_S2_mapping)

MAT_S3_KO_merged <- aggregate(V2 ~ V1, data = MAT_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MAT_S3_mapping$ORF <- rownames(MAT_S3_mapping)
MAT_S3_mapping_KO <- merge(MAT_S3_mapping, MAT_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MAT_S3_mapping_KO) <- MAT_S3_mapping_KO$ORF
MAT_S3_mapping_KO[MAT_S3_mapping_KO == ""] <- NA
rm(MAT_S3_KO, MAT_S3_KO_merged, MAT_S3_mapping)

MAT_S4_KO_merged <- aggregate(V2 ~ V1, data = MAT_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MAT_S4_mapping$ORF <- rownames(MAT_S4_mapping)
MAT_S4_mapping_KO <- merge(MAT_S4_mapping, MAT_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MAT_S4_mapping_KO) <- MAT_S4_mapping_KO$ORF
MAT_S4_mapping_KO[MAT_S4_mapping_KO == ""] <- NA
rm(MAT_S4_KO, MAT_S4_KO_merged, MAT_S4_mapping)

MAT_S5_KO_merged <- aggregate(V2 ~ V1, data = MAT_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MAT_S5_mapping$ORF <- rownames(MAT_S5_mapping)
MAT_S5_mapping_KO <- merge(MAT_S5_mapping, MAT_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MAT_S5_mapping_KO) <- MAT_S5_mapping_KO$ORF
MAT_S5_mapping_KO[MAT_S5_mapping_KO == ""] <- NA
rm(MAT_S5_KO, MAT_S5_KO_merged, MAT_S5_mapping)



MGA_S2_KO_merged <- aggregate(V2 ~ V1, data = MGA_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MGA_S2_mapping$ORF <- rownames(MGA_S2_mapping)
MGA_S2_mapping_KO <- merge(MGA_S2_mapping, MGA_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MGA_S2_mapping_KO) <- MGA_S2_mapping_KO$ORF
MGA_S2_mapping_KO[MGA_S2_mapping_KO == ""] <- NA
rm(MGA_S2_KO, MGA_S2_KO_merged, MGA_S2_mapping)

MGA_S3_KO_merged <- aggregate(V2 ~ V1, data = MGA_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MGA_S3_mapping$ORF <- rownames(MGA_S3_mapping)
MGA_S3_mapping_KO <- merge(MGA_S3_mapping, MGA_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MGA_S3_mapping_KO) <- MGA_S3_mapping_KO$ORF
MGA_S3_mapping_KO[MGA_S3_mapping_KO == ""] <- NA
rm(MGA_S3_KO, MGA_S3_KO_merged, MGA_S3_mapping)

MGA_S4_KO_merged <- aggregate(V2 ~ V1, data = MGA_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MGA_S4_mapping$ORF <- rownames(MGA_S4_mapping)
MGA_S4_mapping_KO <- merge(MGA_S4_mapping, MGA_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MGA_S4_mapping_KO) <- MGA_S4_mapping_KO$ORF
MGA_S4_mapping_KO[MGA_S4_mapping_KO == ""] <- NA
rm(MGA_S4_KO, MGA_S4_KO_merged, MGA_S4_mapping)

MGA_S5_KO_merged <- aggregate(V2 ~ V1, data = MGA_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MGA_S5_mapping$ORF <- rownames(MGA_S5_mapping)
MGA_S5_mapping_KO <- merge(MGA_S5_mapping, MGA_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MGA_S5_mapping_KO) <- MGA_S5_mapping_KO$ORF
MGA_S5_mapping_KO[MGA_S5_mapping_KO == ""] <- NA
rm(MGA_S5_KO, MGA_S5_KO_merged, MGA_S5_mapping)



MOL_S2_KO_merged <- aggregate(V2 ~ V1, data = MOL_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MOL_S2_mapping$ORF <- rownames(MOL_S2_mapping)
MOL_S2_mapping_KO <- merge(MOL_S2_mapping, MOL_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MOL_S2_mapping_KO) <- MOL_S2_mapping_KO$ORF
MOL_S2_mapping_KO[MOL_S2_mapping_KO == ""] <- NA
rm(MOL_S2_KO, MOL_S2_KO_merged, MOL_S2_mapping)

MOL_S3_KO_merged <- aggregate(V2 ~ V1, data = MOL_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MOL_S3_mapping$ORF <- rownames(MOL_S3_mapping)
MOL_S3_mapping_KO <- merge(MOL_S3_mapping, MOL_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MOL_S3_mapping_KO) <- MOL_S3_mapping_KO$ORF
MOL_S3_mapping_KO[MOL_S3_mapping_KO == ""] <- NA
rm(MOL_S3_KO, MOL_S3_KO_merged, MOL_S3_mapping)

MOL_S4_KO_merged <- aggregate(V2 ~ V1, data = MOL_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MOL_S4_mapping$ORF <- rownames(MOL_S4_mapping)
MOL_S4_mapping_KO <- merge(MOL_S4_mapping, MOL_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MOL_S4_mapping_KO) <- MOL_S4_mapping_KO$ORF
MOL_S4_mapping_KO[MOL_S4_mapping_KO == ""] <- NA
rm(MOL_S4_KO, MOL_S4_KO_merged, MOL_S4_mapping)

MOL_S5_KO_merged <- aggregate(V2 ~ V1, data = MOL_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
MOL_S5_mapping$ORF <- rownames(MOL_S5_mapping)
MOL_S5_mapping_KO <- merge(MOL_S5_mapping, MOL_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(MOL_S5_mapping_KO) <- MOL_S5_mapping_KO$ORF
MOL_S5_mapping_KO[MOL_S5_mapping_KO == ""] <- NA
rm(MOL_S5_KO, MOL_S5_KO_merged, MOL_S5_mapping)



PLA_S2_KO_merged <- aggregate(V2 ~ V1, data = PLA_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PLA_S2_mapping$ORF <- rownames(PLA_S2_mapping)
PLA_S2_mapping_KO <- merge(PLA_S2_mapping, PLA_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PLA_S2_mapping_KO) <- PLA_S2_mapping_KO$ORF
PLA_S2_mapping_KO[PLA_S2_mapping_KO == ""] <- NA
rm(PLA_S2_KO, PLA_S2_KO_merged, PLA_S2_mapping)

PLA_S3_KO_merged <- aggregate(V2 ~ V1, data = PLA_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PLA_S3_mapping$ORF <- rownames(PLA_S3_mapping)
PLA_S3_mapping_KO <- merge(PLA_S3_mapping, PLA_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PLA_S3_mapping_KO) <- PLA_S3_mapping_KO$ORF
PLA_S3_mapping_KO[PLA_S3_mapping_KO == ""] <- NA
rm(PLA_S3_KO, PLA_S3_KO_merged, PLA_S3_mapping)

PLA_S4_KO_merged <- aggregate(V2 ~ V1, data = PLA_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PLA_S4_mapping$ORF <- rownames(PLA_S4_mapping)
PLA_S4_mapping_KO <- merge(PLA_S4_mapping, PLA_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PLA_S4_mapping_KO) <- PLA_S4_mapping_KO$ORF
PLA_S4_mapping_KO[PLA_S4_mapping_KO == ""] <- NA
rm(PLA_S4_KO, PLA_S4_KO_merged, PLA_S4_mapping)

PLA_S5_KO_merged <- aggregate(V2 ~ V1, data = PLA_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PLA_S5_mapping$ORF <- rownames(PLA_S5_mapping)
PLA_S5_mapping_KO <- merge(PLA_S5_mapping, PLA_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PLA_S5_mapping_KO) <- PLA_S5_mapping_KO$ORF
PLA_S5_mapping_KO[PLA_S5_mapping_KO == ""] <- NA
rm(PLA_S5_KO, PLA_S5_KO_merged, PLA_S5_mapping)



PON_S2_KO_merged <- aggregate(V2 ~ V1, data = PON_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PON_S2_mapping$ORF <- rownames(PON_S2_mapping)
PON_S2_mapping_KO <- merge(PON_S2_mapping, PON_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PON_S2_mapping_KO) <- PON_S2_mapping_KO$ORF
PON_S2_mapping_KO[PON_S2_mapping_KO == ""] <- NA
rm(PON_S2_KO, PON_S2_KO_merged, PON_S2_mapping)

PON_S3_KO_merged <- aggregate(V2 ~ V1, data = PON_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PON_S3_mapping$ORF <- rownames(PON_S3_mapping)
PON_S3_mapping_KO <- merge(PON_S3_mapping, PON_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PON_S3_mapping_KO) <- PON_S3_mapping_KO$ORF
PON_S3_mapping_KO[PON_S3_mapping_KO == ""] <- NA
rm(PON_S3_KO, PON_S3_KO_merged, PON_S3_mapping)

PON_S4_KO_merged <- aggregate(V2 ~ V1, data = PON_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PON_S4_mapping$ORF <- rownames(PON_S4_mapping)
PON_S4_mapping_KO <- merge(PON_S4_mapping, PON_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PON_S4_mapping_KO) <- PON_S4_mapping_KO$ORF
PON_S4_mapping_KO[PON_S4_mapping_KO == ""] <- NA
rm(PON_S4_KO, PON_S4_KO_merged, PON_S4_mapping)

PON_S5_KO_merged <- aggregate(V2 ~ V1, data = PON_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
PON_S5_mapping$ORF <- rownames(PON_S5_mapping)
PON_S5_mapping_KO <- merge(PON_S5_mapping, PON_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(PON_S5_mapping_KO) <- PON_S5_mapping_KO$ORF
PON_S5_mapping_KO[PON_S5_mapping_KO == ""] <- NA
rm(PON_S5_KO, PON_S5_KO_merged, PON_S5_mapping)



SAU_S2_KO_merged <- aggregate(V2 ~ V1, data = SAU_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
SAU_S2_mapping$ORF <- rownames(SAU_S2_mapping)
SAU_S2_mapping_KO <- merge(SAU_S2_mapping, SAU_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(SAU_S2_mapping_KO) <- SAU_S2_mapping_KO$ORF
SAU_S2_mapping_KO[SAU_S2_mapping_KO == ""] <- NA
rm(SAU_S2_KO, SAU_S2_KO_merged, SAU_S2_mapping)

SAU_S3_KO_merged <- aggregate(V2 ~ V1, data = SAU_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
SAU_S3_mapping$ORF <- rownames(SAU_S3_mapping)
SAU_S3_mapping_KO <- merge(SAU_S3_mapping, SAU_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(SAU_S3_mapping_KO) <- SAU_S3_mapping_KO$ORF
SAU_S3_mapping_KO[SAU_S3_mapping_KO == ""] <- NA
rm(SAU_S3_KO, SAU_S3_KO_merged, SAU_S3_mapping)

SAU_S4_KO_merged <- aggregate(V2 ~ V1, data = SAU_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
SAU_S4_mapping$ORF <- rownames(SAU_S4_mapping)
SAU_S4_mapping_KO <- merge(SAU_S4_mapping, SAU_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(SAU_S4_mapping_KO) <- SAU_S4_mapping_KO$ORF
SAU_S4_mapping_KO[SAU_S4_mapping_KO == ""] <- NA
rm(SAU_S4_KO, SAU_S4_KO_merged, SAU_S4_mapping)

SAU_S5_KO_merged <- aggregate(V2 ~ V1, data = SAU_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
SAU_S5_mapping$ORF <- rownames(SAU_S5_mapping)
SAU_S5_mapping_KO <- merge(SAU_S5_mapping, SAU_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(SAU_S5_mapping_KO) <- SAU_S5_mapping_KO$ORF
SAU_S5_mapping_KO[SAU_S5_mapping_KO == ""] <- NA
rm(SAU_S5_KO, SAU_S5_KO_merged, SAU_S5_mapping)



TON_S2_KO_merged <- aggregate(V2 ~ V1, data = TON_S2_KO, FUN = function(x) paste(unique(x), collapse = "/"))
TON_S2_mapping$ORF <- rownames(TON_S2_mapping)
TON_S2_mapping_KO <- merge(TON_S2_mapping, TON_S2_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(TON_S2_mapping_KO) <- TON_S2_mapping_KO$ORF
TON_S2_mapping_KO[TON_S2_mapping_KO == ""] <- NA
rm(TON_S2_KO, TON_S2_KO_merged, TON_S2_mapping)

TON_S3_KO_merged <- aggregate(V2 ~ V1, data = TON_S3_KO, FUN = function(x) paste(unique(x), collapse = "/"))
TON_S3_mapping$ORF <- rownames(TON_S3_mapping)
TON_S3_mapping_KO <- merge(TON_S3_mapping, TON_S3_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(TON_S3_mapping_KO) <- TON_S3_mapping_KO$ORF
TON_S3_mapping_KO[TON_S3_mapping_KO == ""] <- NA
rm(TON_S3_KO, TON_S3_KO_merged, TON_S3_mapping)

TON_S4_KO_merged <- aggregate(V2 ~ V1, data = TON_S4_KO, FUN = function(x) paste(unique(x), collapse = "/"))
TON_S4_mapping$ORF <- rownames(TON_S4_mapping)
TON_S4_mapping_KO <- merge(TON_S4_mapping, TON_S4_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(TON_S4_mapping_KO) <- TON_S4_mapping_KO$ORF
TON_S4_mapping_KO[TON_S4_mapping_KO == ""] <- NA
rm(TON_S4_KO, TON_S4_KO_merged, TON_S4_mapping)

TON_S5_KO_merged <- aggregate(V2 ~ V1, data = TON_S5_KO, FUN = function(x) paste(unique(x), collapse = "/"))
TON_S5_mapping$ORF <- rownames(TON_S5_mapping)
TON_S5_mapping_KO <- merge(TON_S5_mapping, TON_S5_KO_merged, by.x = "ORF", by.y = "V1", all.x = TRUE)
rownames(TON_S5_mapping_KO) <- TON_S5_mapping_KO$ORF
TON_S5_mapping_KO[TON_S5_mapping_KO == ""] <- NA
rm(TON_S5_KO, TON_S5_KO_merged, TON_S5_mapping)

## Reads mapped against each ORF

# Colonne du nombre de reads mappés sur l'ORF
BRE_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(BRE_S2_mapping_KO[, 8:9]) 
BRE_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(BRE_S3_mapping_KO[, 8:9]) 
BRE_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(BRE_S4_mapping_KO[, 8:9]) 
BRE_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(BRE_S5_mapping_KO[, 8:9]) 

GAB_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(GAB_S2_mapping_KO[, 8:9]) 
GAB_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(GAB_S3_mapping_KO[, 8:9]) 
GAB_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(GAB_S4_mapping_KO[, 8:9]) 
GAB_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(GAB_S5_mapping_KO[, 8:9]) 

MAT_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(MAT_S2_mapping_KO[, 8:9]) 
MAT_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(MAT_S3_mapping_KO[, 8:9]) 
MAT_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(MAT_S4_mapping_KO[, 8:9]) 
MAT_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(MAT_S5_mapping_KO[, 8:9]) 

MGA_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(MGA_S2_mapping_KO[, 8:9]) 
MGA_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(MGA_S3_mapping_KO[, 8:9]) 
MGA_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(MGA_S4_mapping_KO[, 8:9]) 
MGA_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(MGA_S5_mapping_KO[, 8:9]) 

MOL_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(MOL_S2_mapping_KO[, 8:9]) 
MOL_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(MOL_S3_mapping_KO[, 8:9]) 
MOL_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(MOL_S4_mapping_KO[, 8:9]) 
MOL_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(MOL_S5_mapping_KO[, 8:9]) 

PLA_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(PLA_S2_mapping_KO[, 8:9]) 
PLA_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(PLA_S3_mapping_KO[, 8:9]) 
PLA_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(PLA_S4_mapping_KO[, 8:9]) 
PLA_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(PLA_S5_mapping_KO[, 8:9]) 

PON_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(PON_S2_mapping_KO[, 8:9]) 
PON_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(PON_S3_mapping_KO[, 8:9]) 
PON_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(PON_S4_mapping_KO[, 8:9]) 
PON_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(PON_S5_mapping_KO[, 8:9]) 

SAU_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(SAU_S2_mapping_KO[, 8:9]) 
SAU_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(SAU_S3_mapping_KO[, 8:9]) 
SAU_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(SAU_S4_mapping_KO[, 8:9]) 
SAU_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(SAU_S5_mapping_KO[, 8:9]) 

TON_S2_mapping_KO$Nb_Reads_Mapped <- rowSums(TON_S2_mapping_KO[, 8:9]) 
TON_S3_mapping_KO$Nb_Reads_Mapped <- rowSums(TON_S3_mapping_KO[, 8:9]) 
TON_S4_mapping_KO$Nb_Reads_Mapped <- rowSums(TON_S4_mapping_KO[, 8:9]) 
TON_S5_mapping_KO$Nb_Reads_Mapped <- rowSums(TON_S5_mapping_KO[, 8:9]) 

## TPM Calculation
BRE_S2_mapping_KO <- BRE_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_BRE_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_BRE_S2 = sum(RPK_BRE_S2) / 1e6,
         TPM_BRE_S2 = RPK_BRE_S2 / RPK_Sum_per_Million_BRE_S2)
BRE_S3_mapping_KO <- BRE_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_BRE_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_BRE_S3 = sum(RPK_BRE_S3) / 1e6,
         TPM_BRE_S3 = RPK_BRE_S3 / RPK_Sum_per_Million_BRE_S3)
BRE_S4_mapping_KO <- BRE_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_BRE_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_BRE_S4 = sum(RPK_BRE_S4) / 1e6,
         TPM_BRE_S4 = RPK_BRE_S4 / RPK_Sum_per_Million_BRE_S4)
BRE_S5_mapping_KO <- BRE_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_BRE_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_BRE_S5 = sum(RPK_BRE_S5) / 1e6,
         TPM_BRE_S5 = RPK_BRE_S5 / RPK_Sum_per_Million_BRE_S5)

GAB_S2_mapping_KO <- GAB_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_GAB_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_GAB_S2 = sum(RPK_GAB_S2) / 1e6,
         TPM_GAB_S2 = RPK_GAB_S2 / RPK_Sum_per_Million_GAB_S2)
GAB_S3_mapping_KO <- GAB_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_GAB_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_GAB_S3 = sum(RPK_GAB_S3) / 1e6,
         TPM_GAB_S3 = RPK_GAB_S3 / RPK_Sum_per_Million_GAB_S3)
GAB_S4_mapping_KO <- GAB_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_GAB_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_GAB_S4 = sum(RPK_GAB_S4) / 1e6,
         TPM_GAB_S4 = RPK_GAB_S4 / RPK_Sum_per_Million_GAB_S4)
GAB_S5_mapping_KO <- GAB_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_GAB_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_GAB_S5 = sum(RPK_GAB_S5) / 1e6,
         TPM_GAB_S5 = RPK_GAB_S5 / RPK_Sum_per_Million_GAB_S5)


MAT_S2_mapping_KO <- MAT_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MAT_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MAT_S2 = sum(RPK_MAT_S2) / 1e6,
         TPM_MAT_S2 = RPK_MAT_S2 / RPK_Sum_per_Million_MAT_S2)
MAT_S3_mapping_KO <- MAT_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MAT_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MAT_S3 = sum(RPK_MAT_S3) / 1e6,
         TPM_MAT_S3 = RPK_MAT_S3 / RPK_Sum_per_Million_MAT_S3)
MAT_S4_mapping_KO <- MAT_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MAT_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MAT_S4 = sum(RPK_MAT_S4) / 1e6,
         TPM_MAT_S4 = RPK_MAT_S4 / RPK_Sum_per_Million_MAT_S4)
MAT_S5_mapping_KO <- MAT_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MAT_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MAT_S5 = sum(RPK_MAT_S5) / 1e6,
         TPM_MAT_S5 = RPK_MAT_S5 / RPK_Sum_per_Million_MAT_S5)


MGA_S2_mapping_KO <- MGA_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MGA_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MGA_S2 = sum(RPK_MGA_S2) / 1e6,
         TPM_MGA_S2 = RPK_MGA_S2 / RPK_Sum_per_Million_MGA_S2)
MGA_S3_mapping_KO <- MGA_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MGA_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MGA_S3 = sum(RPK_MGA_S3) / 1e6,
         TPM_MGA_S3 = RPK_MGA_S3 / RPK_Sum_per_Million_MGA_S3)
MGA_S4_mapping_KO <- MGA_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MGA_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MGA_S4 = sum(RPK_MGA_S4) / 1e6,
         TPM_MGA_S4 = RPK_MGA_S4 / RPK_Sum_per_Million_MGA_S4)
MGA_S5_mapping_KO <- MGA_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MGA_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MGA_S5 = sum(RPK_MGA_S5) / 1e6,
         TPM_MGA_S5 = RPK_MGA_S5 / RPK_Sum_per_Million_MGA_S5)


MOL_S2_mapping_KO <- MOL_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MOL_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MOL_S2 = sum(RPK_MOL_S2) / 1e6,
         TPM_MOL_S2 = RPK_MOL_S2 / RPK_Sum_per_Million_MOL_S2)
MOL_S3_mapping_KO <- MOL_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MOL_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MOL_S3 = sum(RPK_MOL_S3) / 1e6,
         TPM_MOL_S3 = RPK_MOL_S3 / RPK_Sum_per_Million_MOL_S3)
MOL_S4_mapping_KO <- MOL_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MOL_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MOL_S4 = sum(RPK_MOL_S4) / 1e6,
         TPM_MOL_S4 = RPK_MOL_S4 / RPK_Sum_per_Million_MOL_S4)
MOL_S5_mapping_KO <- MOL_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_MOL_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_MOL_S5 = sum(RPK_MOL_S5) / 1e6,
         TPM_MOL_S5 = RPK_MOL_S5 / RPK_Sum_per_Million_MOL_S5)


PLA_S2_mapping_KO <- PLA_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PLA_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PLA_S2 = sum(RPK_PLA_S2) / 1e6,
         TPM_PLA_S2 = RPK_PLA_S2 / RPK_Sum_per_Million_PLA_S2)
PLA_S3_mapping_KO <- PLA_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PLA_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PLA_S3 = sum(RPK_PLA_S3) / 1e6,
         TPM_PLA_S3 = RPK_PLA_S3 / RPK_Sum_per_Million_PLA_S3)
PLA_S4_mapping_KO <- PLA_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PLA_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PLA_S4 = sum(RPK_PLA_S4) / 1e6,
         TPM_PLA_S4 = RPK_PLA_S4 / RPK_Sum_per_Million_PLA_S4)
PLA_S5_mapping_KO <- PLA_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PLA_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PLA_S5 = sum(RPK_PLA_S5) / 1e6,
         TPM_PLA_S5 = RPK_PLA_S5 / RPK_Sum_per_Million_PLA_S5)


PON_S2_mapping_KO <- PON_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PON_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PON_S2 = sum(RPK_PON_S2) / 1e6,
         TPM_PON_S2 = RPK_PON_S2 / RPK_Sum_per_Million_PON_S2)
PON_S3_mapping_KO <- PON_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PON_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PON_S3 = sum(RPK_PON_S3) / 1e6,
         TPM_PON_S3 = RPK_PON_S3 / RPK_Sum_per_Million_PON_S3)
PON_S4_mapping_KO <- PON_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PON_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PON_S4 = sum(RPK_PON_S4) / 1e6,
         TPM_PON_S4 = RPK_PON_S4 / RPK_Sum_per_Million_PON_S4)
PON_S5_mapping_KO <- PON_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_PON_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_PON_S5 = sum(RPK_PON_S5) / 1e6,
         TPM_PON_S5 = RPK_PON_S5 / RPK_Sum_per_Million_PON_S5)


SAU_S2_mapping_KO <- SAU_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_SAU_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_SAU_S2 = sum(RPK_SAU_S2) / 1e6,
         TPM_SAU_S2 = RPK_SAU_S2 / RPK_Sum_per_Million_SAU_S2)
SAU_S3_mapping_KO <- SAU_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_SAU_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_SAU_S3 = sum(RPK_SAU_S3) / 1e6,
         TPM_SAU_S3 = RPK_SAU_S3 / RPK_Sum_per_Million_SAU_S3)
SAU_S4_mapping_KO <- SAU_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_SAU_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_SAU_S4 = sum(RPK_SAU_S4) / 1e6,
         TPM_SAU_S4 = RPK_SAU_S4 / RPK_Sum_per_Million_SAU_S4)
SAU_S5_mapping_KO <- SAU_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_SAU_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_SAU_S5 = sum(RPK_SAU_S5) / 1e6,
         TPM_SAU_S5 = RPK_SAU_S5 / RPK_Sum_per_Million_SAU_S5)


TON_S2_mapping_KO <- TON_S2_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_TON_S2 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_TON_S2 = sum(RPK_TON_S2) / 1e6,
         TPM_TON_S2 = RPK_TON_S2 / RPK_Sum_per_Million_TON_S2)
TON_S3_mapping_KO <- TON_S3_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_TON_S3 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_TON_S3 = sum(RPK_TON_S3) / 1e6,
         TPM_TON_S3 = RPK_TON_S3 / RPK_Sum_per_Million_TON_S3)
TON_S4_mapping_KO <- TON_S4_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_TON_S4 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_TON_S4 = sum(RPK_TON_S4) / 1e6,
         TPM_TON_S4 = RPK_TON_S4 / RPK_Sum_per_Million_TON_S4)
TON_S5_mapping_KO <- TON_S5_mapping_KO %>%
  mutate(Length_kb = Length / 1000,
         RPK_TON_S5 = Nb_Reads_Mapped / Length_kb,
         RPK_Sum_per_Million_TON_S5 = sum(RPK_TON_S5) / 1e6,
         TPM_TON_S5 = RPK_TON_S5 / RPK_Sum_per_Million_TON_S5)

## Unsure assignation
unsure_BRE_S2 <- BRE_S2_mapping_KO[grepl("/", BRE_S2_mapping_KO$V2), ]
unsure_BRE_S3 <- BRE_S3_mapping_KO[grepl("/", BRE_S3_mapping_KO$V2), ]
unsure_BRE_S4 <- BRE_S4_mapping_KO[grepl("/", BRE_S4_mapping_KO$V2), ]
unsure_BRE_S5 <- BRE_S5_mapping_KO[grepl("/", BRE_S5_mapping_KO$V2), ]

unsure_GAB_S2 <- GAB_S2_mapping_KO[grepl("/", GAB_S2_mapping_KO$V2), ]
unsure_GAB_S3 <- GAB_S3_mapping_KO[grepl("/", GAB_S3_mapping_KO$V2), ]
unsure_GAB_S4 <- GAB_S4_mapping_KO[grepl("/", GAB_S4_mapping_KO$V2), ]
unsure_GAB_S5 <- GAB_S5_mapping_KO[grepl("/", GAB_S5_mapping_KO$V2), ]

unsure_MAT_S2 <- MAT_S2_mapping_KO[grepl("/", MAT_S2_mapping_KO$V2), ]
unsure_MAT_S3 <- MAT_S3_mapping_KO[grepl("/", MAT_S3_mapping_KO$V2), ]
unsure_MAT_S4 <- MAT_S4_mapping_KO[grepl("/", MAT_S4_mapping_KO$V2), ]
unsure_MAT_S5 <- MAT_S5_mapping_KO[grepl("/", MAT_S5_mapping_KO$V2), ]

unsure_MGA_S2 <- MGA_S2_mapping_KO[grepl("/", MGA_S2_mapping_KO$V2), ]
unsure_MGA_S3 <- MGA_S3_mapping_KO[grepl("/", MGA_S3_mapping_KO$V2), ]
unsure_MGA_S4 <- MGA_S4_mapping_KO[grepl("/", MGA_S4_mapping_KO$V2), ]
unsure_MGA_S5 <- MGA_S5_mapping_KO[grepl("/", MGA_S5_mapping_KO$V2), ]

unsure_MOL_S2 <- MOL_S2_mapping_KO[grepl("/", MOL_S2_mapping_KO$V2), ]
unsure_MOL_S3 <- MOL_S3_mapping_KO[grepl("/", MOL_S3_mapping_KO$V2), ]
unsure_MOL_S4 <- MOL_S4_mapping_KO[grepl("/", MOL_S4_mapping_KO$V2), ]
unsure_MOL_S5 <- MOL_S5_mapping_KO[grepl("/", MOL_S5_mapping_KO$V2), ]

unsure_PLA_S2 <- PLA_S2_mapping_KO[grepl("/", PLA_S2_mapping_KO$V2), ]
unsure_PLA_S3 <- PLA_S3_mapping_KO[grepl("/", PLA_S3_mapping_KO$V2), ]
unsure_PLA_S4 <- PLA_S4_mapping_KO[grepl("/", PLA_S4_mapping_KO$V2), ]
unsure_PLA_S5 <- PLA_S5_mapping_KO[grepl("/", PLA_S5_mapping_KO$V2), ]

unsure_PON_S2 <- PON_S2_mapping_KO[grepl("/", PON_S2_mapping_KO$V2), ]
unsure_PON_S3 <- PON_S3_mapping_KO[grepl("/", PON_S3_mapping_KO$V2), ]
unsure_PON_S4 <- PON_S4_mapping_KO[grepl("/", PON_S4_mapping_KO$V2), ]
unsure_PON_S5 <- PON_S5_mapping_KO[grepl("/", PON_S5_mapping_KO$V2), ]

unsure_SAU_S2 <- SAU_S2_mapping_KO[grepl("/", SAU_S2_mapping_KO$V2), ]
unsure_SAU_S3 <- SAU_S3_mapping_KO[grepl("/", SAU_S3_mapping_KO$V2), ]
unsure_SAU_S4 <- SAU_S4_mapping_KO[grepl("/", SAU_S4_mapping_KO$V2), ]
unsure_SAU_S5 <- SAU_S5_mapping_KO[grepl("/", SAU_S5_mapping_KO$V2), ]

unsure_TON_S2 <- TON_S2_mapping_KO[grepl("/", TON_S2_mapping_KO$V2), ]
unsure_TON_S3 <- TON_S3_mapping_KO[grepl("/", TON_S3_mapping_KO$V2), ]
unsure_TON_S4 <- TON_S4_mapping_KO[grepl("/", TON_S4_mapping_KO$V2), ]
unsure_TON_S5 <- TON_S5_mapping_KO[grepl("/", TON_S5_mapping_KO$V2), ]

KO_multiple <- as.data.frame(unique(c(unsure_BRE_S2$V2, unsure_BRE_S3$V2, unsure_BRE_S4$V2, unsure_BRE_S5$V2, 
                                      unsure_GAB_S2$V2, unsure_GAB_S3$V2, unsure_GAB_S4$V2, unsure_GAB_S5$V2,
                                      unsure_MAT_S2$V2, unsure_MAT_S3$V2, unsure_MAT_S4$V2, unsure_MAT_S5$V2,                                            
                                      unsure_MGA_S2$V2, unsure_MGA_S3$V2, unsure_MGA_S4$V2, unsure_MGA_S5$V2, 
                                      unsure_MOL_S2$V2, unsure_MOL_S3$V2, unsure_MOL_S4$V2, unsure_MOL_S5$V2,
                                      unsure_PLA_S2$V2, unsure_PLA_S3$V2, unsure_PLA_S4$V2, unsure_PLA_S5$V2, 
                                      unsure_PON_S2$V2, unsure_PON_S3$V2, unsure_PON_S4$V2, unsure_PON_S5$V2,
                                      unsure_SAU_S2$V2, unsure_SAU_S3$V2, unsure_SAU_S4$V2, unsure_SAU_S5$V2, 
                                      unsure_TON_S2$V2, unsure_TON_S3$V2, unsure_TON_S4$V2, unsure_TON_S5$V2)))

write.csv(KO_multiple, "/Users/floramottet/Documents/Pro/Bio-info/RE_metag/RStudio/KO_multiple.csv")

## KO TPM with ORF merged
BRE_S2_mapping_KO_merged <- BRE_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% # On enlève les colonnes vides
  group_by(V2) %>%
  summarise(TPM_BRE_S2_sum = sum(TPM_BRE_S2))
rownames(BRE_S2_mapping_KO_merged) <- BRE_S2_mapping_KO_merged$V2

BRE_S3_mapping_KO_merged <- BRE_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_BRE_S3_sum = sum(TPM_BRE_S3))
rownames(BRE_S3_mapping_KO_merged) <- BRE_S3_mapping_KO_merged$V2

BRE_S4_mapping_KO_merged <- BRE_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_BRE_S4_sum = sum(TPM_BRE_S4))
rownames(BRE_S4_mapping_KO_merged) <- BRE_S4_mapping_KO_merged$V2

BRE_S5_mapping_KO_merged <- BRE_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_BRE_S5_sum = sum(TPM_BRE_S5))
rownames(BRE_S5_mapping_KO_merged) <- BRE_S5_mapping_KO_merged$V2




GAB_S2_mapping_KO_merged <- GAB_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_GAB_S2_sum = sum(TPM_GAB_S2))
rownames(GAB_S2_mapping_KO_merged) <- GAB_S2_mapping_KO_merged$V2

GAB_S3_mapping_KO_merged <- GAB_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_GAB_S3_sum = sum(TPM_GAB_S3))
rownames(GAB_S3_mapping_KO_merged) <- GAB_S3_mapping_KO_merged$V2

GAB_S4_mapping_KO_merged <- GAB_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_GAB_S4_sum = sum(TPM_GAB_S4))
rownames(GAB_S4_mapping_KO_merged) <- GAB_S4_mapping_KO_merged$V2

GAB_S5_mapping_KO_merged <- GAB_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_GAB_S5_sum = sum(TPM_GAB_S5))
rownames(GAB_S5_mapping_KO_merged) <- GAB_S5_mapping_KO_merged$V2




MAT_S2_mapping_KO_merged <- MAT_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MAT_S2_sum = sum(TPM_MAT_S2))
rownames(MAT_S2_mapping_KO_merged) <- MAT_S2_mapping_KO_merged$V2

MAT_S3_mapping_KO_merged <- MAT_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MAT_S3_sum = sum(TPM_MAT_S3))
rownames(MAT_S3_mapping_KO_merged) <- MAT_S3_mapping_KO_merged$V2

MAT_S4_mapping_KO_merged <- MAT_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MAT_S4_sum = sum(TPM_MAT_S4))
rownames(MAT_S4_mapping_KO_merged) <- MAT_S4_mapping_KO_merged$V2

MAT_S5_mapping_KO_merged <- MAT_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MAT_S5_sum = sum(TPM_MAT_S5))
rownames(MAT_S5_mapping_KO_merged) <- MAT_S5_mapping_KO_merged$V2




MGA_S2_mapping_KO_merged <- MGA_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MGA_S2_sum = sum(TPM_MGA_S2))
rownames(MGA_S2_mapping_KO_merged) <- MGA_S2_mapping_KO_merged$V2

MGA_S3_mapping_KO_merged <- MGA_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MGA_S3_sum = sum(TPM_MGA_S3))
rownames(MGA_S3_mapping_KO_merged) <- MGA_S3_mapping_KO_merged$V2

MGA_S4_mapping_KO_merged <- MGA_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MGA_S4_sum = sum(TPM_MGA_S4))
rownames(MGA_S4_mapping_KO_merged) <- MGA_S4_mapping_KO_merged$V2

MGA_S5_mapping_KO_merged <- MGA_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MGA_S5_sum = sum(TPM_MGA_S5))
rownames(MGA_S5_mapping_KO_merged) <- MGA_S5_mapping_KO_merged$V2




MOL_S2_mapping_KO_merged <- MOL_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MOL_S2_sum = sum(TPM_MOL_S2))
rownames(MOL_S2_mapping_KO_merged) <- MOL_S2_mapping_KO_merged$V2

MOL_S3_mapping_KO_merged <- MOL_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MOL_S3_sum = sum(TPM_MOL_S3))
rownames(MOL_S3_mapping_KO_merged) <- MOL_S3_mapping_KO_merged$V2

MOL_S4_mapping_KO_merged <- MOL_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MOL_S4_sum = sum(TPM_MOL_S4))
rownames(MOL_S4_mapping_KO_merged) <- MOL_S4_mapping_KO_merged$V2

MOL_S5_mapping_KO_merged <- MOL_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_MOL_S5_sum = sum(TPM_MOL_S5))
rownames(MOL_S5_mapping_KO_merged) <- MOL_S5_mapping_KO_merged$V2




PLA_S2_mapping_KO_merged <- PLA_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PLA_S2_sum = sum(TPM_PLA_S2))
rownames(PLA_S2_mapping_KO_merged) <- PLA_S2_mapping_KO_merged$V2

PLA_S3_mapping_KO_merged <- PLA_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PLA_S3_sum = sum(TPM_PLA_S3))
rownames(PLA_S3_mapping_KO_merged) <- PLA_S3_mapping_KO_merged$V2

PLA_S4_mapping_KO_merged <- PLA_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PLA_S4_sum = sum(TPM_PLA_S4))
rownames(PLA_S4_mapping_KO_merged) <- PLA_S4_mapping_KO_merged$V2

PLA_S5_mapping_KO_merged <- PLA_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PLA_S5_sum = sum(TPM_PLA_S5))
rownames(PLA_S5_mapping_KO_merged) <- PLA_S5_mapping_KO_merged$V2




PON_S2_mapping_KO_merged <- PON_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PON_S2_sum = sum(TPM_PON_S2))
rownames(PON_S2_mapping_KO_merged) <- PON_S2_mapping_KO_merged$V2

PON_S3_mapping_KO_merged <- PON_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PON_S3_sum = sum(TPM_PON_S3))
rownames(PON_S3_mapping_KO_merged) <- PON_S3_mapping_KO_merged$V2

PON_S4_mapping_KO_merged <- PON_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PON_S4_sum = sum(TPM_PON_S4))
rownames(PON_S4_mapping_KO_merged) <- PON_S4_mapping_KO_merged$V2

PON_S5_mapping_KO_merged <- PON_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_PON_S5_sum = sum(TPM_PON_S5))
rownames(PON_S5_mapping_KO_merged) <- PON_S5_mapping_KO_merged$V2




SAU_S2_mapping_KO_merged <- SAU_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_SAU_S2_sum = sum(TPM_SAU_S2))
rownames(SAU_S2_mapping_KO_merged) <- SAU_S2_mapping_KO_merged$V2

SAU_S3_mapping_KO_merged <- SAU_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_SAU_S3_sum = sum(TPM_SAU_S3))
rownames(SAU_S3_mapping_KO_merged) <- SAU_S3_mapping_KO_merged$V2

SAU_S4_mapping_KO_merged <- SAU_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_SAU_S4_sum = sum(TPM_SAU_S4))
rownames(SAU_S4_mapping_KO_merged) <- SAU_S4_mapping_KO_merged$V2

SAU_S5_mapping_KO_merged <- SAU_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_SAU_S5_sum = sum(TPM_SAU_S5))
rownames(SAU_S5_mapping_KO_merged) <- SAU_S5_mapping_KO_merged$V2




TON_S2_mapping_KO_merged <- TON_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_TON_S2_sum = sum(TPM_TON_S2))
rownames(TON_S2_mapping_KO_merged) <- TON_S2_mapping_KO_merged$V2

TON_S3_mapping_KO_merged <- TON_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_TON_S3_sum = sum(TPM_TON_S3))
rownames(TON_S3_mapping_KO_merged) <- TON_S3_mapping_KO_merged$V2

TON_S4_mapping_KO_merged <- TON_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_TON_S4_sum = sum(TPM_TON_S4))
rownames(TON_S4_mapping_KO_merged) <- TON_S4_mapping_KO_merged$V2

TON_S5_mapping_KO_merged <- TON_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(TPM_TON_S5_sum = sum(TPM_TON_S5))
rownames(TON_S5_mapping_KO_merged) <- TON_S5_mapping_KO_merged$V2


## KO TPM all samples
samples_list <- list(GAB_S2_mapping_KO_merged, GAB_S3_mapping_KO_merged, GAB_S4_mapping_KO_merged, GAB_S5_mapping_KO_merged, 
                     MAT_S2_mapping_KO_merged, MAT_S3_mapping_KO_merged, MAT_S4_mapping_KO_merged, MAT_S5_mapping_KO_merged, 
                     TON_S2_mapping_KO_merged, TON_S3_mapping_KO_merged, TON_S4_mapping_KO_merged, TON_S5_mapping_KO_merged, 
                     MOL_S2_mapping_KO_merged, MOL_S3_mapping_KO_merged, MOL_S4_mapping_KO_merged, MOL_S5_mapping_KO_merged, 
                     PON_S2_mapping_KO_merged, PON_S3_mapping_KO_merged, PON_S4_mapping_KO_merged, PON_S5_mapping_KO_merged, 
                     MGA_S2_mapping_KO_merged, MGA_S3_mapping_KO_merged, MGA_S4_mapping_KO_merged, MGA_S5_mapping_KO_merged, 
                     PLA_S2_mapping_KO_merged, PLA_S3_mapping_KO_merged, PLA_S4_mapping_KO_merged, PLA_S5_mapping_KO_merged, 
                     BRE_S2_mapping_KO_merged, BRE_S3_mapping_KO_merged, BRE_S4_mapping_KO_merged, BRE_S5_mapping_KO_merged, 
                     SAU_S2_mapping_KO_merged, SAU_S3_mapping_KO_merged, SAU_S4_mapping_KO_merged, SAU_S5_mapping_KO_merged)

KO_Samples_tab <- Reduce(function(x, y) full_join(x, y, by = "V2"),
                         samples_list)

KO_list <- KO_Samples_tab$V2
KO_Samples_tab <- KO_Samples_tab[, -1]
KO_Samples_tab[is.na(KO_Samples_tab)] <- 0
rownames(KO_Samples_tab) <- KO_list

## KO raw counts for DeSeq2 Analysis
BRE_S2_mapping_KO_Raw_Count <- BRE_S2_mapping_KO %>%
filter(!is.na(V2), !grepl("/", V2)) %>% # On enlève les colonnes vides
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_BRE_S2_sum = sum(Nb_Reads_Mapped))
rownames(BRE_S2_mapping_KO_Raw_Count) <- BRE_S2_mapping_KO_Raw_Count$V2

BRE_S3_mapping_KO_Raw_Count <- BRE_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_BRE_S3_sum = sum(Nb_Reads_Mapped))
rownames(BRE_S3_mapping_KO_Raw_Count) <- BRE_S3_mapping_KO_Raw_Count$V2

BRE_S4_mapping_KO_Raw_Count <- BRE_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_BRE_S4_sum = sum(Nb_Reads_Mapped))
rownames(BRE_S4_mapping_KO_Raw_Count) <- BRE_S4_mapping_KO_Raw_Count$V2

BRE_S5_mapping_KO_Raw_Count <- BRE_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_BRE_S5_sum = sum(Nb_Reads_Mapped))
rownames(BRE_S5_mapping_KO_Raw_Count) <- BRE_S5_mapping_KO_Raw_Count$V2




GAB_S2_mapping_KO_Raw_Count <- GAB_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_GAB_S2_sum = sum(Nb_Reads_Mapped))
rownames(GAB_S2_mapping_KO_Raw_Count) <- GAB_S2_mapping_KO_Raw_Count$V2

GAB_S3_mapping_KO_Raw_Count <- GAB_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_GAB_S3_sum = sum(Nb_Reads_Mapped))
rownames(GAB_S3_mapping_KO_Raw_Count) <- GAB_S3_mapping_KO_Raw_Count$V2

GAB_S4_mapping_KO_Raw_Count <- GAB_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_GAB_S4_sum = sum(Nb_Reads_Mapped))
rownames(GAB_S4_mapping_KO_Raw_Count) <- GAB_S4_mapping_KO_Raw_Count$V2

GAB_S5_mapping_KO_Raw_Count <- GAB_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_GAB_S5_sum = sum(Nb_Reads_Mapped))
rownames(GAB_S5_mapping_KO_Raw_Count) <- GAB_S5_mapping_KO_Raw_Count$V2




MAT_S2_mapping_KO_Raw_Count <- MAT_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MAT_S2_sum = sum(Nb_Reads_Mapped))
rownames(MAT_S2_mapping_KO_Raw_Count) <- MAT_S2_mapping_KO_Raw_Count$V2

MAT_S3_mapping_KO_Raw_Count <- MAT_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MAT_S3_sum = sum(Nb_Reads_Mapped))
rownames(MAT_S3_mapping_KO_Raw_Count) <- MAT_S3_mapping_KO_Raw_Count$V2

MAT_S4_mapping_KO_Raw_Count <- MAT_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MAT_S4_sum = sum(Nb_Reads_Mapped))
rownames(MAT_S4_mapping_KO_Raw_Count) <- MAT_S4_mapping_KO_Raw_Count$V2

MAT_S5_mapping_KO_Raw_Count <- MAT_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MAT_S5_sum = sum(Nb_Reads_Mapped))
rownames(MAT_S5_mapping_KO_Raw_Count) <- MAT_S5_mapping_KO_Raw_Count$V2




MGA_S2_mapping_KO_Raw_Count <- MGA_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MGA_S2_sum = sum(Nb_Reads_Mapped))
rownames(MGA_S2_mapping_KO_Raw_Count) <- MGA_S2_mapping_KO_Raw_Count$V2

MGA_S3_mapping_KO_Raw_Count <- MGA_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MGA_S3_sum = sum(Nb_Reads_Mapped))
rownames(MGA_S3_mapping_KO_Raw_Count) <- MGA_S3_mapping_KO_Raw_Count$V2

MGA_S4_mapping_KO_Raw_Count <- MGA_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MGA_S4_sum = sum(Nb_Reads_Mapped))
rownames(MGA_S4_mapping_KO_Raw_Count) <- MGA_S4_mapping_KO_Raw_Count$V2

MGA_S5_mapping_KO_Raw_Count <- MGA_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MGA_S5_sum = sum(Nb_Reads_Mapped))
rownames(MGA_S5_mapping_KO_Raw_Count) <- MGA_S5_mapping_KO_Raw_Count$V2



MOL_S2_mapping_KO_Raw_Count <- MOL_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MOL_S2_sum = sum(Nb_Reads_Mapped))
rownames(MOL_S2_mapping_KO_Raw_Count) <- MOL_S2_mapping_KO_Raw_Count$V2

MOL_S3_mapping_KO_Raw_Count <- MOL_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MOL_S3_sum = sum(Nb_Reads_Mapped))
rownames(MOL_S3_mapping_KO_Raw_Count) <- MOL_S3_mapping_KO_Raw_Count$V2

MOL_S4_mapping_KO_Raw_Count <- MOL_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MOL_S4_sum = sum(Nb_Reads_Mapped))
rownames(MOL_S4_mapping_KO_Raw_Count) <- MOL_S4_mapping_KO_Raw_Count$V2

MOL_S5_mapping_KO_Raw_Count <- MOL_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_MOL_S5_sum = sum(Nb_Reads_Mapped))
rownames(MOL_S5_mapping_KO_Raw_Count) <- MOL_S5_mapping_KO_Raw_Count$V2




PLA_S2_mapping_KO_Raw_Count <- PLA_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PLA_S2_sum = sum(Nb_Reads_Mapped))
rownames(PLA_S2_mapping_KO_Raw_Count) <- PLA_S2_mapping_KO_Raw_Count$V2

PLA_S3_mapping_KO_Raw_Count <- PLA_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PLA_S3_sum = sum(Nb_Reads_Mapped))
rownames(PLA_S3_mapping_KO_Raw_Count) <- PLA_S3_mapping_KO_Raw_Count$V2

PLA_S4_mapping_KO_Raw_Count <- PLA_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PLA_S4_sum = sum(Nb_Reads_Mapped))
rownames(PLA_S4_mapping_KO_Raw_Count) <- PLA_S4_mapping_KO_Raw_Count$V2

PLA_S5_mapping_KO_Raw_Count <- PLA_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PLA_S5_sum = sum(Nb_Reads_Mapped))
rownames(PLA_S5_mapping_KO_Raw_Count) <- PLA_S5_mapping_KO_Raw_Count$V2




PON_S2_mapping_KO_Raw_Count <- PON_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PON_S2_sum = sum(Nb_Reads_Mapped))
rownames(PON_S2_mapping_KO_Raw_Count) <- PON_S2_mapping_KO_Raw_Count$V2

PON_S3_mapping_KO_Raw_Count <- PON_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PON_S3_sum = sum(Nb_Reads_Mapped))
rownames(PON_S3_mapping_KO_Raw_Count) <- PON_S3_mapping_KO_Raw_Count$V2

PON_S4_mapping_KO_Raw_Count <- PON_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PON_S4_sum = sum(Nb_Reads_Mapped))
rownames(PON_S4_mapping_KO_Raw_Count) <- PON_S4_mapping_KO_Raw_Count$V2

PON_S5_mapping_KO_Raw_Count <- PON_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_PON_S5_sum = sum(Nb_Reads_Mapped))
rownames(PON_S5_mapping_KO_Raw_Count) <- PON_S5_mapping_KO_Raw_Count$V2




SAU_S2_mapping_KO_Raw_Count <- SAU_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_SAU_S2_sum = sum(Nb_Reads_Mapped))
rownames(SAU_S2_mapping_KO_Raw_Count) <- SAU_S2_mapping_KO_Raw_Count$V2

SAU_S3_mapping_KO_Raw_Count <- SAU_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_SAU_S3_sum = sum(Nb_Reads_Mapped))
rownames(SAU_S3_mapping_KO_Raw_Count) <- SAU_S3_mapping_KO_Raw_Count$V2

SAU_S4_mapping_KO_Raw_Count <- SAU_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_SAU_S4_sum = sum(Nb_Reads_Mapped))
rownames(SAU_S4_mapping_KO_Raw_Count) <- SAU_S4_mapping_KO_Raw_Count$V2

SAU_S5_mapping_KO_Raw_Count <- SAU_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_SAU_S5_sum = sum(Nb_Reads_Mapped))
rownames(SAU_S5_mapping_KO_Raw_Count) <- SAU_S5_mapping_KO_Raw_Count$V2




TON_S2_mapping_KO_Raw_Count <- TON_S2_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>% 
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_TON_S2_sum = sum(Nb_Reads_Mapped))
rownames(TON_S2_mapping_KO_Raw_Count) <- TON_S2_mapping_KO_Raw_Count$V2

TON_S3_mapping_KO_Raw_Count <- TON_S3_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_TON_S3_sum = sum(Nb_Reads_Mapped))
rownames(TON_S3_mapping_KO_Raw_Count) <- TON_S3_mapping_KO_Raw_Count$V2

TON_S4_mapping_KO_Raw_Count <- TON_S4_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_TON_S4_sum = sum(Nb_Reads_Mapped))
rownames(TON_S4_mapping_KO_Raw_Count) <- TON_S4_mapping_KO_Raw_Count$V2

TON_S5_mapping_KO_Raw_Count <- TON_S5_mapping_KO %>%
  filter(!is.na(V2), !grepl("/", V2)) %>%
  group_by(V2) %>%
  summarise(Nb_Reads_Mapped_TON_S5_sum = sum(Nb_Reads_Mapped))
rownames(TON_S5_mapping_KO_Raw_Count) <- TON_S5_mapping_KO_Raw_Count$V2

samples_list_raw_count <- list(GAB_S2_mapping_KO_Raw_Count, GAB_S3_mapping_KO_Raw_Count, GAB_S4_mapping_KO_Raw_Count, GAB_S5_mapping_KO_Raw_Count, 
                               MAT_S2_mapping_KO_Raw_Count, MAT_S3_mapping_KO_Raw_Count, MAT_S4_mapping_KO_Raw_Count, MAT_S5_mapping_KO_Raw_Count, 
                               TON_S2_mapping_KO_Raw_Count, TON_S3_mapping_KO_Raw_Count, TON_S4_mapping_KO_Raw_Count, TON_S5_mapping_KO_Raw_Count, 
                               MOL_S2_mapping_KO_Raw_Count, MOL_S3_mapping_KO_Raw_Count, MOL_S4_mapping_KO_Raw_Count, MOL_S5_mapping_KO_Raw_Count, 
                               PON_S2_mapping_KO_Raw_Count, PON_S3_mapping_KO_Raw_Count, PON_S4_mapping_KO_Raw_Count, PON_S5_mapping_KO_Raw_Count, 
                               MGA_S2_mapping_KO_Raw_Count, MGA_S3_mapping_KO_Raw_Count, MGA_S4_mapping_KO_Raw_Count, MGA_S5_mapping_KO_Raw_Count, 
                               PLA_S2_mapping_KO_Raw_Count, PLA_S3_mapping_KO_Raw_Count, PLA_S4_mapping_KO_Raw_Count, PLA_S5_mapping_KO_Raw_Count, 
                               BRE_S2_mapping_KO_Raw_Count, BRE_S3_mapping_KO_Raw_Count, BRE_S4_mapping_KO_Raw_Count, BRE_S5_mapping_KO_Raw_Count, 
                               SAU_S2_mapping_KO_Raw_Count, SAU_S3_mapping_KO_Raw_Count, SAU_S4_mapping_KO_Raw_Count, SAU_S5_mapping_KO_Raw_Count)

KO_Samples_tab_raw_count <- Reduce(function(x, y) full_join(x, y, by = "V2"),
                                   samples_list_raw_count)

KO_list_raw_count <- KO_Samples_tab_raw_count$V2
KO_Samples_tab_raw_count <- KO_Samples_tab_raw_count[, -1]
KO_Samples_tab_raw_count[is.na(KO_Samples_tab_raw_count)] <- 0
rownames(KO_Samples_tab_raw_count) <- KO_list_raw_count
colnames(KO_Samples_tab_raw_count) <- Metadata_metaG$Samples

### Récupérer les nombres de reads par voies métaboliques
KO_column_Samples_tab_raw_count <- KO_Samples_tab_raw_count %>%
  rownames_to_column("KO")

KO_Samples_Carbon_Cycle_raw_count <- KO_column_Samples_tab_raw_count[which(KO_column_Samples_tab_raw_count$KO %in% Pathway_KO_KEGG_list$KO),] # On garde les KO qui sont impliqués dans les cycles biogéochimiques d'après DiTing
KO_C_Cycle_list <- KO_Samples_Carbon_Cycle_raw_count$KO
KO_Samples_Carbon_Cycle_raw_count <- KO_Samples_Carbon_Cycle_raw_count[, -1]
rownames(KO_Samples_Carbon_Cycle_raw_count) <- KO_C_Cycle_list



## Save dataset
save.image("/Users/floramottet/Documents/Pro/Bio-info/Data/All_Data_Offiial.RData")