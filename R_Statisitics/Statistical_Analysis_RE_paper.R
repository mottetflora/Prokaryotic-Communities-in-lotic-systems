######################################################################################################
########################### Taxonomic and Functional Statistical Analyses  ###########################
######################################################################################################

##################################################
############### Environment set up ###############
##################################################

#########################
####### Load data #######
#########################
## From Load_Data.R Script
load("/Users/floramottet/Documents/Pro/Bio-info/Data/All_Data_Offiial.RData")

########################
####### Packages #######
########################
library(Hmisc)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(tidyr)
library(vegan)
library(usedist)
library(ape)
library(ggplot2)
library(ggrepel)
library(ade4)
library(adegraphics)
library(tibble)
library(ggVennDiagram)
library(glmmTMB)
library(emmeans)
library(DESeq2)

##########################
###### Color Palets ######
##########################
colors_sidearm <- c("GAB" = "#009E73", "MAT" = "#009E73", "TON" = "#009E73", "MOL" = "#F0E442", "PON" = "#F0E442", "MGA" = "#D55E00", "PLA" = "#D55E00", "BRE" = "#000000", "SAU" = "#000000")
colors_site <- c("GAB-A" = "#009E73", "GAB-B" = "#009E73", "GAB-C" = "#009E73", 
                 "MAT-A" = "#009E73", "MAT-B" = "#009E73", "MAT-C" = "#009E73", 
                 "TON-A" = "#009E73", "TON-B" = "#009E73", "TON-C" = "#009E73", 
                 "MOL-A" = "#F0E442", "MOL-B" = "#F0E442", "MOL-C" = "#F0E442",
                 "PON-A" = "#F0E442", "PON-B" = "#F0E442", "PON-C" = "#F0E442", 
                 "MGA-A" = "#D55E00", 
                 "PLA" = "#D55E00", "BRE" = "#000000", "SAU" = "#000000")
colors_typo <- c("Eupotamon" = "#009E73", "Parapotamon" = "#F0E442", "Plesiopotamon" =  "#D55E00", "Paleopotamon" = "#000000")
colors_campaign <- c("S2" = "#543005", "S3" = "#8C510A", "S4" = "#BF812D", "S5" = "#DFC27D")
colors_typo_boxplot <- c("Eupotamon" = "#009E73", "Parapotamon" = "#F0E442", "Plesiopotamon" =  "#D55E00", "Paleopotamon" = "#252525")

color_typo_season <- c("S2_Eupo" = "#A1D99B","S3_Eupo" = "#74C476", "S4_Eupo" = "#009E73", "S5_Eupo" = "#00441B",
                       "S2_Para" = "#F0E442", "S3_Para" = "#E6AB02", "S4_Para" = "#A6761D", "S5_Para" = "#543005", 
                       "S2_Plesio" = "#FC8D59", "S3_Plesio" = "#D55E00", "S4_Plesio" = "#D7301F", "S5_Plesio" = "#8B0000", 
                       "S2_Paleo" = "#BDBDBD", "S3_Paleo" = "#737373" , "S4_Paleo" = "#525252", "S5_Paleo" = "#000000")

colors_sample <- c("GAB_S2" = "#009E73", "GAB_S3" = "#009E73", "GAB_S4" = "#009E73", "GAB_S5" = "#009E73", 
                   "MAT_S2" = "#009E73", "MAT_S3" = "#009E73", "MAT_S4" = "#009E73", "MAT_S5" = "#009E73", 
                   "TON_S2" = "#009E73", "TON_S3" = "#009E73", "TON_S4" = "#009E73", "TON_S5" = "#009E73", 
                   "MOL_S2" = "#F0E442", "MOL_S3" = "#F0E442", "MOL_S4" = "#F0E442", "MOL_S5" = "#F0E442",
                   "PON_S2" = "#F0E442", "PON_S3" = "#F0E442", "PON_S4" = "#F0E442", "PON_S5" = "#F0E442", 
                   "MGA_S2" = "#D55E00", "MGA_S3" = "#D55E00", "MGA_S4" = "#D55E00", "MGA_S5" = "#D55E00", 
                   "PLA_S2" = "#D55E00", "PLA_S3" = "#D55E00", "PLA_S4" = "#D55E00", "PLA_S5" = "#D55E00", 
                   "BRE_S2" = "#000000", "BRE_S3" = "#000000", "BRE_S4" = "#000000", "BRE_S5" = "#000000",
                   "SAU_S2" = "#000000", "SAU_S3" = "#000000", "SAU_S4" = "#000000", "SAU_S5" = "#000000")

##################################################
############### Taxonomic Analyses ###############
##################################################

#### H1: Hydrological connectivity influence environmental characteristics of sidearms ####
### Environmental variables correlation
## Select only informative data
tab_var_env_temp <- Metadata_used_agreg[,c(1:10, 13:16, 18:19, 21)]
# NA Filtration
tab_var_env_sans_NA <- na.omit(tab_var_env_temp) # 7 sites lost (S3_MAT_B, S4_PLA_A, S4_PLAC, S4_PON_A, S4_PON_B, S4_SAU_C, S5_SAU_C)
Site_Saison_Sans_NA <- tab_var_env_sans_NA$Site_saison
rownames(tab_var_env_sans_NA) <- tab_var_env_sans_NA$Site_saison
tab_var_env <- tab_var_env_sans_NA[,]
rownames(tab_var_env) <- Site_Saison_Sans_NA
tab_var_env$Typologie <- factor(tab_var_env$Typologie, 
                                levels = c("Eupotamon", "Parapotamon", "Plesiopotamon", "Paleopotamon"))
Metadata_used_agreg_pca <- tab_var_env
colnames(Metadata_used_agreg_pca) <- c("Site_season", "O2", "T", "EC", "pH", "NH4", "PO4", "NO3", "P_bioavail", "Clay", "Chl_a", "DOC", "Ignition", "Connectivity", "Campaign", "Site", "Lone")
res_pca_all <- PCA(Metadata_used_agreg_pca %>% select(2:13), scale.unit = TRUE, graph = FALSE)

# Environmental variables contribution
res_pca_all$eig
fviz_pca_var(res_pca_all,
             col.var = "contrib",
             gradient.cols = c("#BDBDBD", "#969696", "#252525"),
             repel = TRUE) +
  scale_x_reverse()
var_pca <- res_pca_all$eig
pc1_var <- round(var_pca[1,2], 1)
pc2_var <- round(var_pca[2,2], 1)   

# Sites representation
coord_sites <- as.data.frame(res_pca_all$ind$coord)
coord_sites$Lone   <- Metadata_used_agreg_pca$Lone
coord_sites$Campaign <- Metadata_used_agreg_pca$Campaign
coord_sites$Connectivity   <- Metadata_used_agreg_pca$Connectivity
coord_sites$Month_typo <- Metadata_used_agreg_pca$Month_Typo
coord_sites$Samples <- paste(coord_sites$Lone, coord_sites$Campaign, sep = "_")
ggplot(coord_sites, aes(x = -Dim.1, y = Dim.2, color = Connectivity)) +
  geom_point(data = coord_sites, aes(color = Connectivity, shape = Campaign), size = 2.5) +
  scale_shape_manual(values = c("S2" = 15, "S3" = 16, "S4" = 17, "S5" = 18),
                     labels = c("May" , "June", "August", "September")) +
  scale_color_manual(values = colors_typo) +
  labs(x = paste0("Axis 1 (", pc1_var, "%)"), 
       y = paste0("Axis 2 (", pc2_var, "%)"), title = "PCA Envionmental Variables") +
  theme_minimal()




