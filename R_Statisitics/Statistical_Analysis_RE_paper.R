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
load("/Users/floramottet/Documents/Pro/Bio-info/Data/All_Data_Official.RData")

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
library(ggplot2)
library(ggrepel)
library(ade4)
library(adegraphics)
library(tibble)
library(ggVennDiagram)
library(doSNOW)
library(itertools)
library(combinat)
library(glmmTMB)
library(emmeans)
library(DESeq2)
library(pheatmap)


##########################
###### Color Palets ######
##########################
colors_typo <- c("Connected" = "#009E73", "Partially connected" = "#F0E442", "Occasionally connected" =  "#D55E00", "Disconnected" = "#000000")

color_typo_season <- c("S2_co" = "#A1D99B","S3_co" = "#74C476", "S4_co" = "#009E73", "S5_co" = "#00441B",
                       "S2_pco" = "#F0E442", "S3_pco" = "#E6AB02", "S4_pco" = "#A6761D", "S5_pco" = "#543005", 
                       "S2_oco" = "#FC8D59", "S3_oco" = "#D55E00", "S4_oco" = "#D7301F", "S5_oco" = "#8B0000", 
                       "S2_dco" = "#BDBDBD", "S3_dco" = "#737373" , "S4_dco" = "#525252", "S5_dco" = "#000000")



colors_sidearm_co_deco <- c("GAB" = "#009E73", "MAT" = "#009E73", "TON" = "#009E73", "MGA" = "#D55E00", "PLA" = "#D55E00", "BRE" = "#D55E00", "SAU" = "#D55E00")


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
###### Functions for taxonomic analysis ######
f <- function(x) log(x+min(x[!is.na(x) & x>1e-12]))

StatBary <-  ggproto("StatBary", Stat, 
                     compute_group = function(data, scales) {
                       tab <- data[c("x", "y")]
                       mns <- colMeans(tab)
                       res <- as.data.frame(as.list(mns))
                       colnames(res) <- c("x", "y")
                       res
                     }, required_aes = c("x", "y"))

stat_bary <- function(mapping = NULL, data = NULL, geom = "point",
                      position = "identity", na.rm = FALSE, show.legend = NA,
                      inherit.aes = TRUE, ...) {
  layer(stat = StatBary, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...))
}

hier.part2 <- function (y, df, gof = "adjR2", parallel = FALSE, 
                        ncores = truc(parallel::detectCores(logical = FALSE)/2),
                        barplot = TRUE){
  pcan <- dim(df)[2]
  if (pcan > 15)
    stop("Number of variables must be < 15 for current implementation",
         call. = FALSE)
  else {
    gfs <- all.rda(y, df, gof = gof, parallel = parallel, ncores = ncores)
    HP <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
    if (barplot) {
      ymin <- min(c(0, floor(min(HP$I.perc) * 0.1) * 10))
      ymax <- ceiling(max(HP$I.perc) * 0.1) * 10
      barplot(t(HP$I.perc), col = c(1), ylim = c(ymin,
                                                 ymax), ylab = "% Independent effects (%I)")
    }
    list(gfs = gfs, IJ = HP$IJ, I.perc = HP$I.perc)
  }
}

all.rda <- function(y, df, gof = "adjR2", condition = NULL, 
                    parallel = TRUE, 
                    ncores = trunc(parallel::detectCores(logical = FALSE)/2)){
  if (!inherits(y, "dist"))
    stop("\ny must be a dist object")
  if (!is.data.frame(df)) 
    stop("df must be a data.frame")
  if (!is.null(condition)) {
    if (!is.character(condition))
      stop("condition must be a character vector")
    if (!condition %in% colnames(df))
      stop("variabile of condition is not in df")
    pcan <- ncol(df)-length(condition)
  }else{
    pcan <- ncol(df)
  }
  n <- (2^pcan) - 1
  combs <- combos1(pcan)$ragged
  if (!is.null(condition)) {
    # place les colonnes des variables de condition à la fin
    df <- cbind(df[, !colnames(df) %in% condition, drop = FALSE], 
                df[, condition, drop = FALSE])
    # rajoute à combs les variable de conditions pour qu'elle soit
    # sytématiquement prise en compte dans le model
    nc <- ncol(df)-length(condition)+1
    mat <- matrix(rep(nc, each = nrow(combs)), nrow(combs))
    combs <- cbind(combs, mat)
  }
  if (!gof %in% c("R2", "adjR2")) {
    cat("\n gof (goodness of fit measure) must equal")
    cat("\n \"R2\" or \"adjR2\" error")
  }
  gfs_null <- 0  # a réfléchir
  if (parallel) {
    cl <- makeCluster(ncores)
    registerDoSNOW(cl)
    clusterExport(cl, c("y", "df", "gof", "current.model2", "condition"),
                  envir = environment())
    gfs <- foreach(it = isplitRows(combs, chunks = ncores), .inorder = TRUE) %dopar% 
      {
        gfsi <- numeric(n <- nrow(it))
        for (i in 1:n) { # a paralleliser
          # le numero de colonne des variabales env
          current.comb <- as.vector(it[i, ][it[i, ] > 0])
          # combn <- paste(names(df)[current.comb], collapse = "")
          new.line <- current.model2(y, current.comb, df,
                                     condition = condition,
                                     gof = gof)
          gfsi[i] <- new.line
        }
        gfsi
      }
    stopCluster(cl)
    gfs <- unlist(gfs, use.names = FALSE)
    gfs <- c(gfs_null, gfs)
  }
  else{
    gfs <- numeric(n+1)
    gfs[1] <- gfs_null
    for (i in 1:n) { # a paralleliser
      if (i%%500 == 0)
        cat(i, "regressions calculated:", n - i, "to go...\n")
      # le numero de colonne des variabales env
      current.comb <- as.vector(combs[i, ][combs[i, ] > 0])
      vars <- names(df)[current.comb]
      if (any(u <- vars%in%condition)){
        vars[u] <- sprintf("Condition(%s)", vars[u])
      }
      # combn <- paste(vars, collapse = "")
      new.line <- current.model2(y, current.comb, df,
                                 condition = condition,
                                 gof = gof)
      gfs[i+1] <- new.line
    }
  }
  
  # if (print.vars) {
  cat("dbrda done: formatting results\n")
  var.names <- character(n+1)
  var.names[1] <- "Theta"
  for (i in 1:n) {
    current.comb <- as.vector(combs[i, ][combs[i, ] >
                                           0])
    vars <- names(data.frame(df)[current.comb])
    if (any(u <- vars %in% condition))
      vars[u] <- sprintf("Condition(%s)", vars[u])
    combi <- paste(vars, collapse = " ")
    var.names[i+1] <- combi
  }
  gfs <- data.frame(`variable combination` = var.names,
                    gof = gfs)
  # }
  gfs
}

current.model2 <- function (y, current.comb, df, condition, gof = "adjR2", 
                            ...){
  comb.data <- df[, current.comb, drop=FALSE]
  colnames(comb.data) <- colnames(df)[current.comb]
  # data <- cbind(y=y, comb.data)
  ny <- deparse(substitute(y))
  # depv <- names(data)[1]
  # n.comb <- ncol(comb.data)
  # xs <- vector("character", n.comb)
  vars <- colnames(comb.data)
  if (any(u <- vars %in%condition)) {
    vars[u] <- sprintf("Condition(%s)", vars[u])
  }
  xs <- paste(vars, collapse = " + ")
  # for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i],
  #                                          "+", sep = "")
  # xs[n.comb] <- names(comb.data)[n.comb]
  # xss <- paste(xs, collapse = " ", sep = "")
  formu <- formula(paste(ny, " ~ ", xs, sep = ""))
  rda <- vegan::dbrda(formu, comb.data, ...)
  r2s <- vegan::RsquareAdj(rda)
  if (gof == "adjR2"){
    gf <- r2s[[2]]
  }else{
    gf <- r2s[[1]]
  }
  gf
}

combos1 <- function (n) {
  #pour avoir l'ensemble des combinaisons possibles
  require(gtools)
  x <- cbind(combinations(n, 1, 1:n), array(0, dim = c(n, n -
                                                         1)))
  for (i in 2:n) {
    nc <- dim(combinations(n, i, 1:n))[1]
    x <- rbind(x, cbind(combinations(n, i, 1:n), array(0,
                                                       dim = c(nc, n - i))))
  }
  len <- dim(x)[1]
  x.index <- cbind(as.vector(1:len), as.vector(x))
  x.index <- cbind(x.index[, 1][x.index[, 2] > 0], x.index[,
                                                           2][x.index[, 2] > 0])
  x.bin <- array(0, dim = c(len, n))
  x.bin[x.index] <- 1
  list(ragged = x, binary = x.bin)
}

partition <- function(nvar, gofs, perc = TRUE, parallel = TRUE,
                      ncores = trunc(parallel::detectCores(logical = FALSE)/2)) {
  # x : le nombre de variables 
  # gofs les mesures de goodness of fits de tous les modèles incluant le modèle null
  # en provenance de all.model (all.rda)
  ragged <- combos1(nvar)$ragged
  r1 <- t(ragged)
  dm <- dim(r1)
  r1 <- as.character(r1)
  r1[r1=="0"] <- ""
  dim(r1) <- dm
  mods <- character(ncol(r1))
  for (i in 1:ncol(r1))
    mods[i] <- paste(r1[, i], collapse = "_", recycle0 = TRUE)
  mods <- sub("_+$", "", mods)
  gofs_null <- gofs[1]
  gofs_wnull <- gofs[-1] 
  x <- 1:nvar
  
  if (parallel){
    ncores <- pmin(ncores, nvar)
    cl <- makeCluster(ncores)
    registerDoSNOW(cl)
    clusterExport(cl, c("mods", "x", "gofs_null", "gofs_wnull"),
                  envir = environment())
    res <- foreach(i = isplitVector(x, chunkSize = 1), .inorder = TRUE,
                   .packages = "combinat") %dopar% {
                     mod1 <- match(i, mods)
                     d_mod1 <- gofs_wnull[mod1]-gofs_null
                     #l'ensemble des chemins possibles = l'ensemble des permutations possibles
                     perms <- permn(setdiff(x, i))
                     N <- length(perms) # le nombre de chemins pour une variable (pourrait être calculer avant)
                     n <- length(perms[[1]])
                     seqs <- 1:n
                     # a chaque chemin il faut récupérer les x models
                     dis <- numeric(N)
                     u <- 1
                     for (j in perms) {
                       mod_2 <- sapply(seqs, function(x) j[1:x])
                       mod_i2 <- mapply(c, i, mod_2)
                       mod_2 <- sapply(mod_2, function(x) paste(sort(x), collapse = "_"))
                       mod_i2 <- sapply(mod_i2, function(x) paste(sort(x), collapse = "_"))
                       
                       gf_mod_wx2 <- gofs_wnull[match(mod_2, mods)]
                       gf_mod2 <- gofs_wnull[match(mod_i2, mods)]
                       
                       dperm <- sum(gf_mod2-gf_mod_wx2)+d_mod1 # ok
                       dperm <- dperm/(n+1) # hierarchy average
                       dis[u] <- dperm
                       u <- u+1
                     }
                     Ic <- mean(dis) # variable average # Ic
                     Jc <- d_mod1 - Ic
                     return(list(Ic = Ic, Jc = Jc))
                   }
    stopCluster(cl)
    Ic <- lapply(res, "[[", "Ic")
    Ic <- unlist(Ic)
    Jc <- lapply(res, "[[", "Jc")
    Jc <- unlist(Jc)
  }else {
    Ic <- Jc <- numeric(nvar)
    v <- 1
    for (i in x) { # i est le numéro de la variable
      # ne pas oublier la première comparaison
      mod1 <- match(i, mods)
      d_mod1 <- gofs_wnull[mod1]-gofs_null
      #l'ensemble des chemins possibles = l'ensemble des permutations possibles
      perms <- permn(setdiff(x, i))
      N <- length(perms) # le nombre de chemins pour une variable (pourrait être calculer avant)
      n <- length(perms[[1]])
      seqs <- 1:n
      # a chaque chemin il faut récupérer les x models
      dis <- numeric(N)
      u <- 1
      for (j in perms) {
        mod_2 <- sapply(seqs, function(x) j[1:x])
        mod_i2 <- mapply(c, i, mod_2)
        mod_2 <- sapply(mod_2, function(x) paste(sort(x), collapse = "_"))
        mod_i2 <- sapply(mod_i2, function(x) paste(sort(x), collapse = "_"))
        
        gf_mod_wx2 <- gofs_wnull[match(mod_2, mods)]
        gf_mod2 <- gofs_wnull[match(mod_i2, mods)]
        
        dperm <- sum(gf_mod2-gf_mod_wx2)+d_mod1 # ok
        dperm <- dperm/(n+1) # hierarchy average
        dis[u] <- dperm
        u <- u+1
      }
      Ic[v] <- mean(dis) # variable average # Ic
      Jc[v] <- d_mod1 - Ic[v]
      v <- v+1
    }
  }
  if (perc) {
    Ic <- Ic/sum(Ic)
    Jc <- Jc/sum(Jc)
  }
  return(list(I = Ic, J = Jc))
}

###########

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
                                levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))
Metadata_used_agreg_pca <- tab_var_env
colnames(Metadata_used_agreg_pca) <- c("Site_season", "DO", "T", "EC", "pH", "NH4", "PO4", "NO3", "Pbiodispo", "Clay", "Chlî", "DOC", "OM", "Connectivity", "Campaign", "Site", "Lone")
# Transforming asymetric variables to reduce extremes values influences
Metadata_used_agreg_pca[, c("Chlî", "DOC", "EC", "OM",
                 "NH4", "NO3", "Pbiodispo", "PO4")] <- 
  lapply(Metadata_used_agreg_pca[, c("Chlî", "DOC", "EC", "OM",
                          "NH4", "NO3", "Pbiodispo", "PO4")], f)
res_pca_all <- PCA(Metadata_used_agreg_pca[, c(2:13)], scale.unit = TRUE, graph = FALSE)

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
###########

#### H2: Hydrological connectivity influences prokaryotic communities structure ####
Table_ab_sum_rel <- decostand(Table_ab_sum, method = "total", MARGIN = 1)
Table_ab_sum_rel <-
  cbind(Table_ab_sum_rel,
        Metadata[match(rownames(Table_ab_sum_rel), sub("-1$", "", Metadata$Pcr)),
                 c("Lone", "Nom", "Typologie", "Saison", "Site")])
## Bray-Curtis matrix
bray_dist_stand_data_ab_sum <- vegdist(Table_ab_sum_rel[,!names(Table_ab_sum_rel) %in% c("Lone","Saison", "Typologie", "Site", "Nom")], method="bray")
## PCoA on mean communities by site
pcoa_com_moy <- dudi.pco(bray_dist_stand_data_ab_sum, scannf = FALSE, nf = 3)
pts_com_moy <- as.data.frame(pcoa_com_moy$li[,1:2])
colnames(pts_com_moy) <- c("PCoA1", "PCoA2")
eig_pcoa <- pcoa_com_moy$eig
eig_pcoa_pct <- eig_pcoa / sum(eig_pcoa[eig_pcoa > 0]) * 100
xlab_pcoa <- paste0("PCoA1 (", round(eig_pcoa_pct[1], 1), "%)")
ylab_pcoa <- paste0("PCoA2 (", round(eig_pcoa_pct[2], 1), "%)")
pts_com_moy <-
  cbind(pts_com_moy,
        Metadata_used[match(rownames(pts_com_moy), sub("-1$", "", Metadata$Pcr)),
                 c("Lone", "Nom", "Typologie", "Saison", "Site")])

ggplot(pts_com_moy, aes(x = PCoA1, y = PCoA2, color = Typologie)) +
  geom_point(data = pts_com_moy, aes(color = Typologie, shape = Saison), size = 2.5) +
  stat_ellipse(aes(color = Typologie, group = Typologie), level = 0.95) +
  scale_shape_manual(values = c("S2" = 15, "S3" = 16, "S4" = 17, "S5" = 18), labels = c("May" , "June", "August", "September")) +
  scale_color_manual(values = colors_typo) +
  coord_cartesian(xlim = c(-0.65, 0.65), ylim = c(-0.4, 0.4)) +
  labs(x = xlab_pcoa, y = ylab_pcoa, title = "PCoA - Bray Curtis ~ Typologie + Mois") +
  theme_minimal()

## perMANOVA
# Centroids calculs to avoid pseudo-replication
dist_centroids_all <- dist_multi_centroids(bray_dist_stand_data_ab_sum,
                                           Metadata$Site_saison, squared = FALSE)
# Hydrological connectivity test
print("#PERMANOVA to test variations between hydrological connectivities")
perma_typo_com_moy <- adonis2(dist_centroids_all ~ Typologie, data = Metadata_used_agreg, by = "margin")
print(perma_typo_com_moy)
###########

#### H3: Hydrological connectivity promotes anaerobic taxa ####
## Specific richness
Lone_site <- Metadata$Site
richesse_spe_site <- matrix(NA, nrow = 291, ncol = 1000)
rownames(richesse_spe_site) <- Metadata$Site_saison
richesse_moyenne_site <- matrix(NA, nrow = 25, ncol = 1000)
rownames(richesse_moyenne_site) <- unique(Metadata$Site)
richesse_sd_site <- matrix(NA, nrow = 25, ncol = 1000)
rownames(richesse_sd_site) <- unique(Metadata$Site)
# 1000 Rarefactions to have a mean and no influence of the rarefaction : Can take 10 min
for(i in 1:1000){
  # Rarefaction
  Data_raref_spec_rich <- rrarefy(round(Table_ab_sum, 0), min(round(rowSums(Table_ab_sum), 0)))
  # Part of specific richness
  richesse_spe_site[, i] <- specnumber(Data_raref_spec_rich)
  richesse_moyenne_site[, i] <- tapply(richesse_spe_site[, i], Lone_site, mean)
  richesse_sd_site[, i] <- tapply(richesse_spe_site[, i], Lone_site, sd)
}
table_richesse_moy_site <- data.frame(
  Richesse_moyenne = apply(richesse_moyenne_site, 1, mean, na.rm = TRUE),
  SD = apply(richesse_sd_site, 1, sd, na.rm = TRUE)
)
table_richesse_moy_site$Site <- rownames(table_richesse_moy_site)
table_richesse_moy_site$Site <- factor(table_richesse_moy_site$Site, 
                                       levels = c("GAB-A", "GAB-B", "GAB-C", "MAT-A", "MAT-B", "MAT-C","TON-A", "TON-B", "TON-C",  
                                                  "MOL-A", "MOL-B", "MOL-C", "PON-A", "PON-B", "PON-C", 
                                                  "MGA-A", "PLA-A", "PLA-B", "PLA-C", 
                                                  "BRE-A", "BRE-B", "BRE-C", "SAU-A", "SAU-B", "SAU-C"))
table_richesse_moy_site$Lone <- substr(table_richesse_moy_site$Site, 1, 3)
table_richesse_moy_site$Lone <- factor(table_richesse_moy_site$Lone, 
                                       levels = c("GAB", "MAT", "TON", "MOL", "PON", "MGA", "PLA", "BRE", "SAU"))
table_richesse_moy_site$Typo <- c(rep("Disconnected", 3), rep("Connected", 5), "Occasionally connected", 
                                  rep("Partially connected", 3), rep("Occasionally connected", 3), 
                                  rep("Partially connected", 3), rep("Disconnected", 3), rep("Connected", 4))
table_richesse_moy_site$Typo <- factor(table_richesse_moy_site$Typo, 
                                       levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))

# Generalized Linear Mixed Model
print("# GLMM and Chi2 test for specific richness ~ Typologie")
spec_rich_model <- glmmTMB(Richesse_moyenne ~ Typo + (1 | Lone) + (1 | Site), 
                           family = poisson(link = "log"),
                           data = table_richesse_moy_site
)
spec_rich_null_model <- glmmTMB(Richesse_moyenne ~ 1 + (1 | Lone) + (1 | Site),
                                family = poisson(link="log"),
                                data = table_richesse_moy_site
)
anova(spec_rich_null_model, spec_rich_model)
emmeans(spec_rich_model, pairwise ~ Typo)
table_richesse_moy_site %>%
  group_by(Typo) %>%
  summarise(
    moyenne = mean(Richesse_moyenne),
    ecart_type = sd(Richesse_moyenne),
    n = n()
  )

## Shannon and Simpson indices
Lone_site <- Metadata$Site
shannon_site <- matrix(NA, nrow = 291, ncol = 1000)
rownames(shannon_site) <- Metadata$Site_saison
shannon_moy_site <- matrix(NA, nrow = 25, ncol = 1000)
rownames(shannon_moy_site) <- unique(Metadata$Site)
shannon_sd_site <- matrix(NA, nrow = 25, ncol = 1000)
rownames(shannon_sd_site) <- unique(Metadata$Site)
simpson_site <- matrix(NA, nrow = 291, ncol = 1000)
rownames(simpson_site) <- Metadata$Site_saison
simpson_moy_site <- matrix(NA, nrow = 25, ncol = 1000)
rownames(simpson_moy_site) <- unique(Metadata$Site)
simpson_sd_site <- matrix(NA, nrow = 25, ncol = 1000)
rownames(simpson_sd_site) <- unique(Metadata$Site)

# 1000 Rarefactions to have a mean and no influence of the rarefaction : Can take 40 min
for(i in 1:1000){
  # Rarefaction
  Data_raref_div_index <- rrarefy(round(Table_ab_sum, 0), min(round(rowSums(Table_ab_sum), 0)))
  # Part of Shannon index
  shannon_site[, i] <- diversity(Data_raref_div_index, index = "shannon")
  shannon_moy_site[, i] <- tapply(shannon_site[, i], Lone_site, mean)
  shannon_sd_site[, i] <- tapply(shannon_site[, i], Lone_site, sd)
  # Part of Simpson index
  simpson_site[, i] <- diversity(Data_raref_div_index, index = "simpson")
  simpson_moy_site[, i] <- tapply(simpson_site[, i], Lone_site, mean)
  simpson_sd_site[, i] <- tapply(simpson_site[, i], Lone_site, sd)
}
table_shannon_moy <- data.frame(
  Shannon_moyen = apply(shannon_moy_site, 1, mean, na.rm = TRUE),
  SD = apply(shannon_sd_site, 1, sd, na.rm = TRUE)
)
table_simpson_moy <- data.frame(
  Simpson_moy = apply(simpson_moy_site, 1, mean, na.rm = TRUE),
  SD = apply(simpson_sd_site, 1, sd, na.rm = TRUE)
)
table_shannon_moy$Site <- rownames(table_shannon_moy)
table_shannon_moy$Site <- factor(table_shannon_moy$Site, 
                                 levels = c("GAB-A", "GAB-B", "GAB-C", "MAT-A", "MAT-B", "MAT-C","TON-A", "TON-B", "TON-C",  
                                            "MOL-A", "MOL-B", "MOL-C", "PON-A", "PON-B", "PON-C", 
                                            "MGA-A", "PLA-A", "PLA-B", "PLA-C", 
                                            "BRE-A", "BRE-B", "BRE-C", "SAU-A", "SAU-B", "SAU-C"))
table_shannon_moy$Lone <- substr(table_shannon_moy$Site, 1, 3)
table_shannon_moy$Lone <- factor(table_shannon_moy$Lone, 
                                 levels = c("GAB", "MAT", "TON", "MOL", "PON", "MGA", "PLA", "BRE", "SAU"))
table_shannon_moy$Typo <- c(rep("Disconnected", 3), rep("Connected", 5), "Occasionally connected", 
                            rep("Partially connected", 3), rep("Occasionally connected", 3), 
                            rep("Partially connected", 3), rep("Disconnected", 3), rep("Connected", 4))
table_shannon_moy$Typo <- factor(table_shannon_moy$Typo, 
                                 levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))
table_simpson_moy$Site <- rownames(table_simpson_moy)
table_simpson_moy$Site <- factor(table_simpson_moy$Site, 
                                 levels = c("GAB-A", "GAB-B", "GAB-C", "MAT-A", "MAT-B", "MAT-C","TON-A", "TON-B", "TON-C",  
                                            "MOL-A", "MOL-B", "MOL-C", "PON-A", "PON-B", "PON-C", 
                                            "MGA-A", "PLA-A", "PLA-B", "PLA-C", 
                                            "BRE-A", "BRE-B", "BRE-C", "SAU-A", "SAU-B", "SAU-C"))
table_simpson_moy$Lone <- substr(table_simpson_moy$Site, 1, 3)
table_simpson_moy$Lone <- factor(table_simpson_moy$Lone, 
                                 levels = c("GAB", "MAT", "TON", "MOL", "PON", "MGA", "PLA", "BRE", "SAU"))
table_simpson_moy$Typo <- c(rep("Disconnected", 3), rep("Connected", 5), "Occasionally connected", 
                            rep("Partially connected", 3), rep("Occasionally connected", 3), 
                            rep("Partially connected", 3), rep("Disconnected", 3), rep("Connected", 4))
table_simpson_moy$Typo <- factor(table_simpson_moy$Typo, 
                                 levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))

# Generalized Linear Mixed Model
print("# LMM and Chi2 test for shannon index ~ Typologie")
shannon_model <- glmmTMB(Shannon_moyen ~ Typo + (1 | Lone) + (1 | Site),
                         family = gaussian(link = "identity"),
                         data = table_shannon_moy
)
shannon_model_null <- glmmTMB(Shannon_moyen ~ 1 + (1 | Lone) + (1 | Site),
                              family = gaussian(link = "identity"),
                              data = table_shannon_moy
)
anova(shannon_model_null, shannon_model)
emmeans(shannon_model, pairwise ~ Typo)
table_shannon_moy %>%
  group_by(Typo) %>%
  summarise(
    moyenne = mean(Shannon_moyen),
    ecart_type = sd(Shannon_moyen),
    n = n()
  )
print("# LMM and Chi2 test for Simpson index ~ Typologie")
simpson_model <- glmmTMB(Simpson_moy ~ Typo + (1 | Lone) + (1 | Site),
                         family = gaussian(link = "identity"),
                         data = table_simpson_moy
)
simpson_model_null <- glmmTMB(Simpson_moy ~ Typo + (1 | Lone) + (1 | Site), 
                              family = nbinom2(link = "log"),
                              data = table_simpson_moy
)
anova(simpson_model_null, simpson_model)
emmeans(simpson_model, pairwise ~ Typo)
table_simpson_moy %>%
  group_by(Typo) %>%
  summarise(
    moyenne = mean(Simpson_moy),
    ecart_type = sd(Simpson_moy),
    n = n()
  )

## Phyla Proportions
bacteries <- Assigned_OTU[Assigned_OTU$Domaine == "Bacteria", ]
archees <- Assigned_OTU[Assigned_OTU$Domaine == "Archaea", ]
pct_Data_raref_per_site <- Data_raref_per_site * 100 / rowSums(Data_raref_per_site)
table(Assigned_OTU$Domaine)
pct_archees_per_site <- pct_Data_raref_per_site[, colnames(pct_Data_raref_per_site) %in% rownames(archees)]
pct_bact_per_site <- pct_Data_raref_per_site[, colnames(pct_Data_raref_per_site) %in% rownames(bacteries)]
## Histogram per phylum
# Remove unassigned data at domain level (450 OTUs)
Assigned_Domain <- Assigned_OTU[Assigned_OTU$Domaine != "Other",]
Assigned_Domain$Domain_Phylum <- paste(Assigned_Domain$Domaine, Assigned_Domain$Phylum, sep = "_")
pct_phylum_transpo_per_site <- pct_Data_raref_per_site[, colnames(pct_Data_raref_per_site) %in% rownames(Assigned_Domain)]
pct_phylum_per_site <- rowsum(t(pct_phylum_transpo_per_site), group = Assigned_Domain$Domain_Phylum)
sort(rowSums(pct_phylum_per_site)/25, decreasing = TRUE)[1:15] # Up to 1%
phylum_ab_list_per_site <- names(sort(rowSums(pct_phylum_per_site)/25, decreasing = TRUE)[1:15])
pct_phylum_ab_transpo_per_site <- pct_phylum_per_site[rownames(pct_phylum_per_site) %in% phylum_ab_list_per_site,]
pct_phylum_other_per_site <- pct_phylum_per_site[!rownames(pct_phylum_per_site) %in% phylum_ab_list_per_site,]
pct_phylum_ab_transpo_per_site <- rbind(pct_phylum_ab_transpo_per_site, rep(NA, ncol(pct_phylum_ab_transpo_per_site)))
rownames(pct_phylum_ab_transpo_per_site)[16] <- "Archaea_Other"
pct_phylum_ab_transpo_per_site <- rbind(pct_phylum_ab_transpo_per_site, rep(NA, ncol(pct_phylum_ab_transpo_per_site)))
rownames(pct_phylum_ab_transpo_per_site)[17] <- "Bacteria_Other"
pct_phylum_ab_transpo_per_site[16,] <- colSums(pct_phylum_other_per_site[1:10,])
pct_phylum_ab_transpo_per_site[17,] <- colSums(pct_phylum_other_per_site[11:73,])
pct_phylum_ab_per_site <- as.data.frame(t(pct_phylum_ab_transpo_per_site))
pct_phylum_ab_per_site$Typo <- c(rep("Disconnected", 3), 
                                 rep("Connected", 6), "Occasionally connected",
                                 rep("Partially connected", 3),
                                 rep("Occasionally connected", 3),
                                 rep("Partially connected", 3),
                                 rep("Disconnected", 3), 
                                 rep("Connected", 3))
phylum_long_per_site <- pct_phylum_ab_per_site %>%
  rownames_to_column(var = "Site") %>%
  pivot_longer(
    cols = colnames(pct_phylum_ab_per_site[,1:17]),
    names_to = "Phylum",
    values_to = "Proportion"
  )
phylum_long_per_site$Typo <- factor(phylum_long_per_site$Typo,
                                    levels = c("Disconnected", "Occasionally connected", "Partially connected", "Connected"))
phylum_long_per_site$Site <- factor(phylum_long_per_site$Site, 
                                    levels = c("SAU-A", "SAU-B", "SAU-C",
                                               "BRE-A", "BRE-B", "BRE-C",
                                               "PLA-A", "PLA-B", "PLA-C", 
                                               "MGA-A",
                                               "PON-A", "PON-B", "PON-C", 
                                               "MOL-A", "MOL-B", "MOL-C",
                                               "TON-A", "TON-B", "TON-C",
                                               "MAT-A", "MAT-B", "MAT-C",
                                               "GAB-A", "GAB-B", "GAB-C"))
phylum_long_per_site$Phylum <- factor(phylum_long_per_site$Phylum, 
                                      levels = c("Archaea_Thermoplasmatota", "Archaea_Thermoproteota", "Archaea_Other",
                                                 "Bacteria_Acidobacteriota", "Bacteria_Bacteroidota", "Bacteria_Candidatus", "Bacteria_Chloroflexota", "Bacteria_Ignavibacteriota", "Bacteria_Latescibacterota", 
                                                 "Bacteria_Nitrospirota", "Bacteria_Planctomycetota", "Bacteria_Pseudomonadota", "Bacteria_Spirochaetota", 
                                                 "Bacteria_Sva0485", "Bacteria_Thermodesulfobacteriota", "Bacteria_Verrucomicrobiota", "Bacteria_Other"))
ggplot(phylum_long_per_site, aes(x = Site, y = Proportion, fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = c("Archaea_Thermoplasmatota" = "#8E0152", "Archaea_Thermoproteota" = "#DE77AE", "Archaea_Other" = "#FDE0EF",
                               "Bacteria_Acidobacteriota" = "#C6DBEF", "Bacteria_Bacteroidota" ="#4EB3D3", "Bacteria_Candidatus" = "#6BAED6",  "Bacteria_Chloroflexota" = "#4292C6", "Bacteria_Ignavibacteriota" = "#2171B5", "Bacteria_Latescibacterota" = "#0066CC", 
                               "Bacteria_Nitrospirota" = "#08519C", "Bacteria_Planctomycetota" = "#08306B", "Bacteria_Pseudomonadota" = "#00441B", "Bacteria_Spirochaetota" = "#006D2C", 
                               "Bacteria_Sva0485" = "#238B45", "Bacteria_Thermodesulfobacteriota" = "#41AE76", "Bacteria_Verrucomicrobiota" = "#66C2A4", "Bacteria_Other" = "#99D8C9")) +
  labs(x = "Site", y = "Proportion of reads (%)", fill = "Phylum") +
  facet_wrap(~ Typo, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.spacing = unit(1, "lines")
  )
## Proportions of every class
sort(rowSums(rowsum(t(pct_Data_raref_per_site[, colnames(pct_Data_raref_per_site) %in% rownames(Assigned_Domain)]), group = Assigned_Domain$Classe))/25, decreasing = TRUE)
## Proportions of every order
sort(rowSums(rowsum(t(pct_Data_raref_per_site[, colnames(pct_Data_raref_per_site) %in% rownames(Assigned_Domain)]), group = Assigned_Domain$Ordre))/25, decreasing = TRUE)
###########

#### H4: Environmental variables explaining prokaryotic communities structure ####
# Transform asymetric variables in log + min(x) values
Metadata_temp <- Metadata
Metadata_temp[c(29, 31, 16, 32, 19, 21, 22, 20)] <- 
  lapply(Metadata_temp[, c(29, 31, 16, 32, 19, 21, 22, 20)], f)
Metadata_temp <- Metadata_temp[order(Metadata_temp$Site_saison), ]

# Identify which samples contains NA in metadata
samples_na <- sapply(lapply(Metadata_temp[, c(13, 15:17, 19:22, 26, 29, 31, 32)], is.na), which) %>%
  Reduce(union, .) %>% 
  sort()

# Remove samples with NA in metadata
Metadata_temp <- Metadata_temp[-samples_na, ]
tab_ab_sum_rel_db_rda <- Table_ab_sum_rel[!sub("\\d$", "", rownames(Table_ab_sum_rel)) %in%
                                            Metadata_temp$Site_saison[samples_na], ]

# Keep the 5 last columns and remove OTUs with lowest abundances (abundance < e-12)
auxi <- tab_ab_sum_rel_db_rda[, -((ncol(tab_ab_sum_rel_db_rda)-4):ncol(tab_ab_sum_rel_db_rda))]
auxi <- auxi[, colSums(auxi)>1e-12]
tab_ab_sum_rel_db_rda <- cbind(auxi, tab_ab_sum_rel_db_rda[, ((ncol(tab_ab_sum_rel_db_rda)-4):ncol(tab_ab_sum_rel_db_rda))])

# Centroids per site
bray_dist_stand_db_rda <- as.matrix(vegdist(tab_ab_sum_rel_db_rda[,!names(tab_ab_sum_rel_db_rda) %in% c("Lone","Saison", "Typologie", "Site", "Nom")], method="bray"))
Metadata_temp <- Metadata_temp[match(rownames(bray_dist_stand_db_rda), sub("-\\d$", "", Metadata_temp$Pcr)), ]
dist_centroids_all_db_rda <- dist_multi_centroids(bray_dist_stand_db_rda,
                                                  Metadata_temp$Site_saison, squared = FALSE)
Metadata_temp_bis <- Metadata_temp[match(attr(dist_centroids_all_db_rda, "Labels"), 
                                         Metadata_temp$Site_saison),]

# Informative metadata
Metadata_db_rda <- cbind(Metadata_temp_bis[, c(13, 15:17, 19:22, 26, 29, 31, 32)],
             Metadata_temp_bis[, c("Saison"), drop = FALSE])

# Variance Inflation Factor 
vif.cca(capscale(dist_centroids_all_db_rda ~ ., data = Metadata_db_rda[, -13])) 

# partial db-RDA models
rda_all_seasons_minp <- capscale(dist_centroids_all_db_rda ~ Condition(Saison),
                                 data = Metadata_db_rda)
rda_all_seasons_maxp <- 
  capscale(dist_centroids_all_db_rda ~ Condition(Saison) + O2_mg + T + cond + pH + NH4_col + PO4_col + 
             NO3_col + Pbiodispo_bentho + clay_pct + chla_col + COD_col +
             paf_bentho, data = Metadata_db_rda)
step_dbrdap <- ordiR2step(rda_all_seasons_minp,
                          scope = rda_all_seasons_maxp,
                          direction = "both",
                          R2scope = TRUE,
                          pstep = 1000,
                          trace = FALSE)
step_dbrdap
varsp <- attr(terms(step_dbrdap), "term.labels")
varsp <- sub("(^Condition\\(|)", "", varsp)
varsp <- sub("\\)", "", varsp)

# db-RDA formula after variables selction
formu <- sprintf("%s ~ %s",  deparse(substitute(dist_centroids_all_db_rda)),
                 paste(attr(terms(step_dbrdap), "term.labels"), collapse = " + "))
formu <- as.formula(formu)

# Associated statistics
rda_all_seasons_p <- dbrda(formu, Metadata_db_rda)
rda_all_seasons_p
centroid_dist_rda <- capscale(formula = dist_centroids_all_db_rda ~ Condition(Saison) + COD_col + paf_bentho + PO4_col + chla_col + clay_pct + T + O2_mg + NH4_col, data = Metadata_db_rda) 
RsquareAdj(centroid_dist_rda) 
anova(centroid_dist_rda, permutations = 999)
anova.cca(centroid_dist_rda, by = "margin")
anova.cca(centroid_dist_rda, by = "axis")

## Graph
text_dbRDA1_all <- paste("dbRDA1 (",(4.4254/13.1435)*100," %)")
text_dbRDA2_all <- paste("dbRDA2 (",(0.9561/13.1435)*100," %)")
site_all <- as.data.frame(scores(centroid_dist_rda)[2])
variables_all <-as.data.frame(scores(centroid_dist_rda)[4])
rownames(variables_all) <- c("DOC", "OM", "PO4", "Chlî", "Clay", "T", "DO", "NH4")
site_all$Lone <- Metadata_temp_bis$Lone
site_all$Month_Connect <- Metadata_temp_bis$Month_typo
site_all$Typo <- Metadata_temp_bis$Typologie
site_all$Month_Connect <- factor(site_all$Month_Connect, 
                                 levels = c("S2_co","S3_co", "S4_co", "S5_co",
                                            "S2_pco", "S3_pco", "S4_pco", "S5_pco", 
                                            "S2_oco", "S3_oco", "S4_oco", "S5_oco", 
                                            "S2_dco", "S3_dco", "S4_dco", "S5_dco"))
site_all$Lone <- factor(site_all$Lone, 
                        levels = c("BRE", "GAB", "MAT", "MGA", "MOL", "PLA", "PON", "SAU", "TON"))
site_all$Typo <- factor(site_all$Typo, levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))

ggplot(site_all, aes(x = sites.CAP1, y = -sites.CAP2)) +
  geom_point(data=site_all, size = 4, aes(x=sites.CAP1, y=-sites.CAP2, color = Month_Connect, shape = Lone)) +
  labs(col="Hydrological connectivity") +
  geom_segment(data=variables_all,aes(x=biplot.CAP1*0, y=-biplot.CAP2*0, xend=biplot.CAP1*2, yend=-biplot.CAP2*2), arrow = arrow(length = unit(0.15, "cm"))) +
  geom_text_repel(data=variables_all,aes(x=biplot.CAP1*2,y=-biplot.CAP2*2,label=rownames(variables_all))) +
  scale_color_manual(values = color_typo_season) +
  scale_shape_manual(values = c("BRE" = 9, "GAB" = 15, "MAT" = 7,
                                "MGA" = 17, "MOL" = 10, "PLA" = 6, 
                                "PON" = 16, "SAU" = 18, "TON" = 12)) +
  labs(x=text_dbRDA1_all, y=text_dbRDA2_all) +
  ggtitle("db-RDA for All Seasons") +
  theme(plot.title = element_text(size=10)) +
  theme_bw()

# Hierarchical partitioning
all1p <- all.rda(dist_centroids_all_db_rda, Metadata_db_rda[, varsp],
                 gof = "adjR2", condition = "Saison", 
                 parallel = TRUE)
hps6 <- partition(length(varsp)-1, all1p$gof, 
                  parallel = TRUE, ncores = length(varsp)-1)

hier_part_tab <- data.frame(setdiff(varsp, "Saison"), do.call(cbind, hps6))
colnames(hier_part_tab) <- c("Variables", "I", "J")
###########

#### H5: Some taxa are found only in one hydrological connectivity ####
# All OTUs
typo_lone <- data.frame(
  Site = rownames(Data_raref_per_lone),
  Typologie = c(
    "Disconnected", "Connected","Connected", "Occasionally connected", "Partially connected", "Occasionally connected", "Partially connected", "Disconnected", "Connected"
  )
)
Data_raref_co <- colnames(Data_raref_per_lone)[colSums(Data_raref_per_lone[c(2, 3, 9),]) > 0]
Data_raref_pco <- colnames(Data_raref_per_lone)[colSums(Data_raref_per_lone[c(5, 7),]) > 0] 
Data_raref_oco <- colnames(Data_raref_per_lone)[colSums(Data_raref_per_lone[c(4, 6),]) > 0] 
Data_raref_dco <- colnames(Data_raref_per_lone)[colSums(Data_raref_per_lone[c(1, 8),]) > 0]
Data_raref_per_lone_typo <- list("Connected" = Data_raref_co, 
                                 "Partially connected" = Data_raref_pco, 
                                 "Occasionally connected" = Data_raref_oco, 
                                 "Disconnected" = Data_raref_dco)
ggVennDiagram(
  Data_raref_per_lone_typo,
  label_alpha = 0,
  set_color = c("#009E73", "#F0E442", "#D55E00", "#000000")
) +
  scale_fill_gradient(low = "white", high = "grey")
# Number of reads per OTUs
rowSums(Data_raref_per_lone)
colSums(Data_raref_per_lone)
nb_reads_per_sample <- 7930
OTUs_co <- Data_raref_per_lone[colSums(Data_raref_per_lone[c(2, 3, 9),]) > 0]
OTUs_pco <- Data_raref_per_lone[colSums(Data_raref_per_lone[c(5, 7),]) > 0] 
OTUs_oco <- Data_raref_per_lone[colSums(Data_raref_per_lone[c(4, 6),]) > 0] 
OTUs_dco <- Data_raref_per_lone[colSums(Data_raref_per_lone[c(1, 8),]) > 0]
OTUs_tot <- Data_raref_per_lone[colSums(Data_raref_per_lone) > 0] # 22 433 OTUs abundances > 0
OTUs_tot <- Data_raref_per_site_per_season[colSums(Data_raref_per_site_per_season) > 0]
OTUs_all_names <- Reduce(intersect, list(
  colnames(OTUs_co),
  colnames(OTUs_pco),
  colnames(OTUs_oco),
  colnames(OTUs_dco)
))
reads_OTUs_all_names <- rowSums(Data_raref_per_lone[, which(colnames(Data_raref_per_lone) %in% OTUs_all_names)])
sum(reads_OTUs_all_names)
(sum(reads_OTUs_all_names)*100) / (nb_reads_per_sample*9)
# Connected (co)
OTUs_only_co <- Reduce(setdiff, list(
  colnames(OTUs_co),
  c(colnames(OTUs_pco), colnames(OTUs_oco), colnames(OTUs_dco))
))
reads_OTUs_only_co <- rowSums(Data_raref_per_lone[, which(colnames(Data_raref_per_lone) %in% OTUs_only_co)])
sum(reads_OTUs_only_co)
(sum(reads_OTUs_only_co)*100) / (nb_reads_per_sample*9)
(sum(reads_OTUs_only_co)*100) / (nb_reads_per_sample*3)
# Partially connected (pco)
OTUs_only_pco <- Reduce(setdiff, list(
  colnames(OTUs_pco),
  c(colnames(OTUs_co), colnames(OTUs_oco), colnames(OTUs_dco))
))
reads_OTUs_only_pco <- rowSums(Data_raref_per_lone[, which(colnames(Data_raref_per_lone) %in% OTUs_only_pco)])
sum(reads_OTUs_only_pco)
(sum(reads_OTUs_only_pco)*100) / (nb_reads_per_sample*9)
(sum(reads_OTUs_only_pco)*100) / (nb_reads_per_sample*2)
# Occasionally connected (oco)
OTUs_only_oco <- Reduce(setdiff, list(
  colnames(OTUs_oco),
  c(colnames(OTUs_co), colnames(OTUs_pco), colnames(OTUs_paleo))
))
reads_OTUs_only_oco <- rowSums(Data_raref_per_lone[, which(colnames(Data_raref_per_lone) %in% OTUs_only_oco)])
sum(reads_OTUs_only_oco)
(sum(reads_OTUs_only_oco)*100) / (nb_reads_per_sample*9)
(sum(reads_OTUs_only_oco)*100) / (nb_reads_per_sample*2)
# Disconnected (dco)
OTUs_only_dco <- Reduce(setdiff, list(
  colnames(OTUs_dco),
  c(colnames(OTUs_co), colnames(OTUs_pco), colnames(OTUs_oco))
))
reads_OTUs_only_dco <- rowSums(Data_raref_per_lone[, which(colnames(Data_raref_per_lone) %in% OTUs_only_dco)])
sum(reads_OTUs_only_dco)
(sum(reads_OTUs_only_dco)*100) / (nb_reads_per_sample*9)
(sum(reads_OTUs_only_dco)*100) / (nb_reads_per_sample*2)


###########

#################################################
############## Functional Analyses ##############
#################################################
#### H6: Hydrological connectivity influences prokaryotic functional profiles structure ####
## Environmental Data
tab_var_env_temp_all_KO <- Metadata_metaG[,c(10:18, 22:23)]
tab_var_env_all_KO_sans_NA <- na.omit(tab_var_env_temp_all_KO) 
Site_Saison_Sans_NA_all_KO <- rownames(tab_var_env_all_KO_sans_NA)
tab_var_env_all_KO_used <- tab_var_env_all_KO_sans_NA # PON_S4 and PLA_S4 lost (pH missing)
## Experimental design
tab_saison_site <- Metadata_metaG[which(rownames(Metadata_metaG) %in% Site_Saison_Sans_NA_all_KO), c(4:6)]
tab_saison_site$Sampling_campaign <- as.factor(tab_saison_site$Sampling_campaign)
tab_saison_site$Lone <- as.factor(tab_saison_site$Lone)
tab_saison_site$Lone <- factor(tab_saison_site$Lone, 
                               levels = c("GAB", "MAT", "TON", "MOL", "PON", "MGA", "PLA", "BRE", "SAU"))
tab_saison_site$Connectivity <- factor(tab_saison_site$Connectivity, 
                                       levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))
## KO data
tab_compt_temp <- as.data.frame(t(KO_Samples_tab))
rownames(tab_compt_temp) <- substr(rownames(tab_compt_temp), 5, 10)
tab_compt_temp_sans_NA <- tab_compt_temp[which(rownames(tab_compt_temp) %in% Site_Saison_Sans_NA_all_KO),]
## Standardization of Abundances
tab_compt <- decostand(tab_compt_temp_sans_NA, "hellinger")
KO_TPM_per_sample <- t(KO_Samples_tab)
stand_KO_TPM <- decostand(KO_TPM_per_sample, "hellinger")
bray_KO_TPM <- vegdist(stand_KO_TPM, method = "bray")
## PCoA with all sites
pcoa_KO_bray <- dudi.pco(bray_KO_TPM, scannf = FALSE, nf = 3)
pts_KO_bray <- as.data.frame(pcoa_KO_bray$li[,1:2])
colnames(pts_KO_bray) <- c("PCoA1", "PCoA2")
eig_KO_bray <- pcoa_KO_bray$eig
eig_KO_bray_pct <- eig_KO_bray / sum(eig_KO_bray[eig_KO_bray > 0]) * 100
xlab_KO_bray <- paste0("PCoA1 (", round(eig_KO_bray_pct[1], 1), "%)")
ylab_KO_bray <- paste0("PCoA2 (", round(eig_KO_bray_pct[2], 1), "%)")

pts_KO_bray <-
  cbind(pts_KO_bray,
        Metadata_metaG[match(substr(rownames(pts_KO_bray), 5, 10), rownames(Metadata_metaG)),
                      c("Samples", "Name", "Connectivity", "Sampling_campaign")])
rownames(pts_KO_bray) <- pts_KO_bray$Samples

pts_KO_bray$Season_form <- factor(
  ifelse(pts_KO_bray$Sampling_campaign == "S2", "circle", 
         ifelse(pts_KO_bray$Sampling_campaign == "S3", "diamond", 
                ifelse(pts_KO_bray$Sampling_campaign == "S4", "square", "triangle"))))

ggplot(pts_KO_bray, aes(x = -PCoA1, y = PCoA2, color = Connectivity)) +
  geom_point(data = pts_KO_bray, aes(color = Connectivity, shape = Sampling_campaign), size = 2.5) +
  stat_ellipse(aes(color = Connectivity, group = Connectivity), type = "norm", level = 0.67) +
  scale_shape_manual(values = c("S2" = 15, "S3" = 16, "S4" = 17, "S5" = 18), 
                     labels = c("May" , "June", "August", "September")) +
  scale_color_manual(values = colors_typo) +
  coord_cartesian(xlim = c(-0.4, 0.3), ylim = c(-0.15, 0.2)) +
  labs(x = xlab_KO_bray, y = ylab_KO_bray, title = "PCoA - Bray Curtis ~ Typologie + Mois") +
  theme_minimal()

## perMANOVA
# Hydrological connectivity test
print("#PERMANOVA to test variations between hydrological connectivities")
perma_typo_KO <- adonis2(bray_KO_TPM ~ Connectivity, data = Metadata_metaG, by = "margin")
print(perma_typo_KO)

## PCoA for KOs involved in Carbon cycle
tab_compt_Carbon_temp <- tab_compt_temp[, which(colnames(tab_compt_temp) %in% Pathway_KO_KEGG_list$KO)]
tab_compt_C <- decostand(tab_compt_Carbon_temp, "hellinger")
stand_KO_TPM_C <- decostand(tab_compt_C, "hellinger")
bray_KO_TPM_C <- vegdist(stand_KO_TPM_C, method = "bray")
pcoa_KO_bray_C <- dudi.pco(bray_KO_TPM_C, scannf = FALSE, nf = 3)
pts_KO_bray_C <- as.data.frame(pcoa_KO_bray_C$li[,1:2])
colnames(pts_KO_bray_C) <- c("PCoA1", "PCoA2")
eig_KO_bray_C <- pcoa_KO_bray_C$eig
eig_KO_bray_C_pct <- eig_KO_bray_C / sum(eig_KO_bray_C[eig_KO_bray_C > 0]) * 100
xlab_KO_bray_C <- paste0("PCoA1 (", round(eig_KO_bray_C_pct[1], 1), "%)")
ylab_KO_bray_C <- paste0("PCoA2 (", round(eig_KO_bray_C_pct[2], 1), "%)")

pts_KO_bray_C <-
  cbind(pts_KO_bray_C,
        Metadata_metaG[match(rownames(pts_KO_bray_C), rownames(Metadata_metaG)),
                       c("Samples", "Name", "Connectivity", "Sampling_campaign")])
rownames(pts_KO_bray_C) <- pts_KO_bray_C$Samples

pts_KO_bray_C$Season_form <- factor(
  ifelse(pts_KO_bray_C$Sampling_campaign == "S2", "circle", 
         ifelse(pts_KO_bray_C$Sampling_campaign == "S3", "diamond", 
                ifelse(pts_KO_bray_C$Sampling_campaign == "S4", "square", "triangle"))))

ggplot(pts_KO_bray_C, aes(x = -PCoA1, y = PCoA2, color = Connectivity)) +
  geom_point(data = pts_KO_bray_C, aes(color = Connectivity, shape = Sampling_campaign), size = 2.5) +
  stat_ellipse(aes(color = Connectivity, group = Connectivity), type = "norm", level = 0.67) +
  scale_shape_manual(values = c("S2" = 15, "S3" = 16, "S4" = 17, "S5" = 18), 
                     labels = c("May" , "June", "August", "September")) +
  scale_color_manual(values = colors_typo) +
  coord_cartesian(xlim = c(-0.4, 0.3), ylim = c(-0.15, 0.2)) +
  labs(x = xlab_KO_bray_C, y = ylab_KO_bray_C, title = "PCoA - Bray Curtis ~ Typologie + Mois") +
  theme_minimal()

print("#PERMANOVA to test variations between hydrological connectivities")
perma_typo_KO_C_Cycle <- adonis2(bray_KO_TPM_C ~ Connectivity, data = Metadata_metaG, by = "margin")
print(perma_typo_KO_C_Cycle)
###########

#### H7: Environmental variables explaining prokaryotic functional profiles structure ####
# Transform asymetric variables in log + min(x) values
Metadata_metaG_temp <- Metadata_metaG

Metadata_metaG_temp[c(11, 13:16, 20:22)] <- 
  lapply(Metadata_metaG_temp[, c(11, 13:16, 20:22)], f)
Metadata_metaG_temp <- Metadata_metaG_temp[order(Metadata_metaG_temp$Samples), ]

# Which samples have NAs ? 
samples_metaG_na <- sapply(lapply(Metadata_metaG_temp[, c(9:17, 20:22)], is.na), which) %>%
  Reduce(union, .) %>% 
  sort()

# Remove samples with NA in metadata
Metadata_metaG_temp <- Metadata_metaG_temp[-samples_metaG_na, ]
KO_TPM_per_sample_db_rda <- KO_TPM_per_sample[!substr(rownames(KO_TPM_per_sample), 5, 10) %in%
                                            Metadata_metaG_temp$Samples[samples_metaG_na], ]
stand_KO_TPM_db_rda <- decostand(KO_TPM_per_sample_db_rda, "hellinger")
bray_KO_TPM_db_rda <- vegdist(stand_KO_TPM_db_rda, method = "bray")

# Informative metadata
Metadata_metaG_db_rda <- cbind(Metadata_metaG_temp[, c(9:17, 20:22)],
                         Metadata_metaG_temp[, c("Sampling_campaign"), drop = FALSE])

# Variance Inflation Factor 
vif.cca(capscale(bray_KO_TPM_db_rda ~ ., data = Metadata_metaG_db_rda[, -13])) 

# partial db-RDA models
rda_metaG_minp <- capscale(bray_KO_TPM_db_rda ~ Condition(Sampling_campaign),
                                 data = Metadata_metaG_db_rda)
rda_metaG_maxp <- capscale(bray_KO_TPM_db_rda ~ Condition(Sampling_campaign) + O2_mg + T + cond + pH + NH4_col + PO4_col + NO3_col + Pbiodispo_bentho + clay_pct + chla_col + COD_col +paf_bentho, data = Metadata_metaG_db_rda)


step_metaG_dbrdap <- ordiR2step(rda_metaG_minp,
                          scope = rda_metaG_maxp,
                          direction = "both",
                          R2scope = TRUE,
                          pstep = 1000,
                          trace = FALSE)
step_metaG_dbrdap
varsp_metaG <- attr(terms(step_metaG_dbrdap), "term.labels")
varsp_metaG <- sub("(^Condition\\(|)", "", varsp_metaG)
varsp_metaG <- sub("\\)", "", varsp_metaG)

# db-RDA formula after variables selction
formu_metaG <- sprintf("%s ~ %s",  deparse(substitute(bray_KO_TPM_db_rda)),
                 paste(attr(terms(step_metaG_dbrdap), "term.labels"), collapse = " + "))
formu_metaG <- as.formula(formu_metaG)

# Associated statistics
rda_metaG_p <- dbrda(formu_metaG, Metadata_metaG_db_rda)
rda_metaG_p
centroid_dist_rda_metaG <- capscale(formula = bray_KO_TPM_db_rda ~ Condition(Sampling_campaign), data = Metadata_metaG_db_rda) 
RsquareAdj(centroid_dist_rda_metaG) 
anova(centroid_dist_rda_metaG, permutations = 999)
anova.cca(centroid_dist_rda_metaG, by = "margin")
anova.cca(centroid_dist_rda_metaG, by = "axis")

## Graph ==> Modèle avec la rda du metabar
text_dbRDA1_all <- paste("dbRDA1 (",(4.4254/13.1435)*100," %)")
text_dbRDA2_all <- paste("dbRDA2 (",(0.9561/13.1435)*100," %)")
site_all <- as.data.frame(scores(centroid_dist_rda)[2])
variables_all <-as.data.frame(scores(centroid_dist_rda)[4])
rownames(variables_all) <- c("DOC", "OM", "PO4", "Chlî", "Clay", "T", "DO", "NH4")
site_all$Lone <- Metadata_temp_bis$Lone
site_all$Month_Connect <- Metadata_temp_bis$Month_typo
site_all$Typo <- Metadata_temp_bis$Typologie
site_all$Month_Connect <- factor(site_all$Month_Connect, 
                                 levels = c("S2_co","S3_co", "S4_co", "S5_co",
                                            "S2_pco", "S3_pco", "S4_pco", "S5_pco", 
                                            "S2_oco", "S3_oco", "S4_oco", "S5_oco", 
                                            "S2_dco", "S3_dco", "S4_dco", "S5_dco"))
site_all$Lone <- factor(site_all$Lone, 
                        levels = c("BRE", "GAB", "MAT", "MGA", "MOL", "PLA", "PON", "SAU", "TON"))
site_all$Typo <- factor(site_all$Typo, levels = c("Connected", "Partially connected", "Occasionally connected", "Disconnected"))

ggplot(site_all, aes(x = sites.CAP1, y = -sites.CAP2)) +
  geom_point(data=site_all, size = 4, aes(x=sites.CAP1, y=-sites.CAP2, color = Month_Connect, shape = Lone)) +
  labs(col="Hydrological connectivity") +
  geom_segment(data=variables_all,aes(x=biplot.CAP1*0, y=-biplot.CAP2*0, xend=biplot.CAP1*2, yend=-biplot.CAP2*2), arrow = arrow(length = unit(0.15, "cm"))) +
  geom_text_repel(data=variables_all,aes(x=biplot.CAP1*2,y=-biplot.CAP2*2,label=rownames(variables_all))) +
  scale_color_manual(values = color_typo_season) +
  scale_shape_manual(values = c("BRE" = 9, "GAB" = 15, "MAT" = 7,
                                "MGA" = 17, "MOL" = 10, "PLA" = 6, 
                                "PON" = 16, "SAU" = 18, "TON" = 12)) +
  labs(x=text_dbRDA1_all, y=text_dbRDA2_all) +
  ggtitle("db-RDA for All Seasons") +
  theme(plot.title = element_text(size=10)) +
  theme_bw()

# Hierarchical partitioning
all1p <- all.rda(dist_centroids_all_db_rda, Metadata_db_rda[, varsp],
                 gof = "adjR2", condition = "Saison", 
                 parallel = TRUE)
hps6 <- partition(length(varsp)-1, all1p$gof, 
                  parallel = TRUE, ncores = length(varsp)-1)

hier_part_tab <- data.frame(setdiff(varsp, "Saison"), do.call(cbind, hps6))
colnames(hier_part_tab) <- c("Variables", "I", "J")
###########

#### H8: Is there a differential presence of some genes  cycle pathways abundances depend on hydrological connectivity? ####
## Number of reads per pathway
KO_column_Samples_tab_raw_count <- KO_Samples_tab_raw_count %>%
  rownames_to_column("KO")
KO_Samples_Carbon_Cycle_raw_count <- KO_column_Samples_tab_raw_count[which(KO_column_Samples_tab_raw_count$KO %in% Pathway_KO_KEGG_list$KO),]
KO_C_Cycle_list <- KO_Samples_Carbon_Cycle_raw_count$KO
KO_Samples_Carbon_Cycle_raw_count <- KO_Samples_Carbon_Cycle_raw_count[, -1]
KO_Samples_Carbon_Cycle_raw_count_sans_inter <- KO_Samples_Carbon_Cycle_raw_count[,which(colnames(KO_Samples_Carbon_Cycle_raw_count) %in% rownames(Metadata_DESeq))]
KO_Samples_Carbon_Cycle_raw_count_sans_inter <- KO_Samples_Carbon_Cycle_raw_count_sans_inter[, order(colnames(KO_Samples_Carbon_Cycle_raw_count_sans_inter))]
KO_Samples_Carbon_Cycle_raw_count_sans_inter <- as.data.frame(KO_Samples_Carbon_Cycle_raw_count_sans_inter)
rownames(KO_Samples_Carbon_Cycle_raw_count_sans_inter) <- KO_C_Cycle_list

## DESeq2 analysis
deseq_data <- DESeqDataSetFromMatrix(
  countData = KO_Samples_Carbon_Cycle_raw_count_sans_inter,
  colData = Metadata_DESeq,
  design = ~ Typo_court
)
deseq_data_final <- DESeq(deseq_data)
resultsNames(deseq_data_final) # Connected est au numérateur, Disconnected est au dénominateur
deseq_result <- results(deseq_data_final)
deseq_foldchange_signif <- deseq_result[which(deseq_result$padj < 0.05 & abs(deseq_result$log2FoldChange) > 2),]
nrow(deseq_foldchange_signif)
deseq_foldchange_signif[order(abs(deseq_foldchange_signif$log2FoldChange), decreasing = TRUE), ]
deseq_foldchange_signif_connected <- subset(deseq_foldchange_signif, log2FoldChange > 2)
nrow(deseq_foldchange_signif_connected)
deseq_foldchange_signif_disconnected <- subset(deseq_foldchange_signif, log2FoldChange < -2)
nrow(deseq_foldchange_signif_disconnected)
## Graphe
deseq_KO <- matrix(NA, nrow = 45, ncol = 2)
colnames(deseq_KO) <- c("KO", "Carbon_role")
deseq_KO <- as.data.frame(deseq_KO)
deseq_KO$KO <- rownames(deseq_foldchange_signif)
deseq_KO$Carbon_role <- c("Methanogenesis", "Carbon fixation", rep("Methanogenesis", 8), "Respiration", "Respiration", "Methanogenesis", "Respiration", "Central metabolism", "Respiration", "Central metabolism", rep("Methanogenesis", 3), "Methanotrophy", rep("Methanogenesis", 3), "Carbon fixation", "Central metabolism", rep("Methanogenesis", 3), "Central metabolism", rep("Methanogenesis", 5), "Central metabolism", rep("Methanogenesis", 4), "Central metabolism", "Carbon fixation", rep("Methanogenesis", 2), "Central metabolism")
plot(deseq_result$log2FoldChange,
     -log10(deseq_result$pvalue),
     pch = ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Central metabolism" ],
                  17, 
                  ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Respiration" ],
                         15, 
                         ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Carbon fixation" ],
                                18, 
                                ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Methanogenesis" ],
                                       8, 
                                       ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Methanotrophy"],
                                              16, 21))))),
     col = ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Central metabolism" ],
                  "black", 
                  ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Respiration" ],
                         "#084594", 
                         ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Carbon fixation" ],
                                "#6BAED6", 
                                ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Methanogenesis" ],
                                       "#A50F15", 
                                       ifelse(rownames(deseq_result) %in% deseq_KO$KO[deseq_KO$Carbon_role == "Methanotrophy"],
                                              "#EF8A62", "grey"))))),
     cex = 1.5,
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)")
abline(v = c(-2, 2), col = "#404040", lty = 4)
abline(h = -log10(0.05), col = "#404040", lty = 4)
## KO abundances (RPMS)
KO_mean_ab <- aggregate(
  t(KO_Samples_Carbon_Cycle_raw_count_sans_inter),
  by = list(Metadata_DESeq$Lone),
  FUN = mean
)
KO_mean_ab <- KO_mean_ab[, -1]
rownames(KO_mean_ab) <- unique(Metadata_DESeq$Lone)
pheatmap(log10(t(KO_mean_ab[,which(colnames(KO_mean_ab) %in% deseq_KO$KO)]) + 1))

###########

#### H9: Carbon cycle pathways abundances depend on hydrological connectivity ####
KO_column_Samples_tab <- KO_Samples_tab %>%
  rownames_to_column("KO")
KO_Samples_C_Cycle <- KO_column_Samples_tab[which(KO_column_Samples_tab$KO %in% Pathway_KO_KEGG_list$KO),]
KO_Samples_long <- KO_Samples_C_Cycle %>%
  pivot_longer(-KO, names_to = "Sample", values_to = "Count")
KO_multiple_pathways <- Pathway_KO_KEGG_list %>%
  group_by(KO) %>%
  summarise(n_pathways = n())
## Weighting abundances by number of pathways involving each KO
KO_Samples_long_weighted <- KO_Samples_long %>%
  inner_join(Pathway_KO_KEGG_list, by = "KO") %>%
  inner_join(KO_multiple_pathways, by = "KO") %>%
  mutate(weighted_abundance = Count / n_pathways)
## Pathways in carbon cycle
Pathway_abundance_long <- KO_Samples_long_weighted %>%
  group_by(Module, Carbon_cycle_section, Carbon_pathway, Sample) %>%
  summarise(Abundance_pathway = sum(weighted_abundance), .groups = "drop") 
Pathway_abundance <- Pathway_abundance_long[which(Pathway_abundance_long$Carbon_pathway != "Other"),] %>%
  pivot_wider(names_from = Sample, values_from = Abundance_pathway)
Pathway_abundance <- Pathway_abundance[order(match(Pathway_abundance$Carbon_pathway, unique(Pathway_KO_KEGG_list$Carbon_pathway))),] 
Carbon_cycle_abundance <- as.data.frame(t(Pathway_abundance[,-c(1:2)]))
colnames(Carbon_cycle_abundance) <- Carbon_cycle_abundance[1,]
Carbon_cycle_abundance <- as.data.frame(Carbon_cycle_abundance[-1,])
rownames(Carbon_cycle_abundance) <- substr(rownames(Carbon_cycle_abundance), 5, 10)
Carbon_cycle_abundance <- Carbon_cycle_abundance[,]%>%
  mutate(across(everything(), as.numeric))
Carbon_pathways <- colnames(Carbon_cycle_abundance)
long_data <- Carbon_cycle_abundance %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    -Sample,
    names_to = "Pathway",
    values_to = "Abundance"
  )
long_data$Sample <- factor(long_data$Sample, 
                           levels = c("GAB_S2", "GAB_S3", "GAB_S4", "GAB_S5", "MAT_S2", "MAT_S3", "MAT_S4", "MAT_S5", "TON_S2", "TON_S3", "TON_S4", "TON_S5", 
                                      "MOL_S2", "MOL_S3", "MOL_S4", "MOL_S5", "PON_S2", "PON_S3", "PON_S4", "PON_S5", 
                                      "MGA_S2", "MGA_S3", "MGA_S4", "MGA_S5", "PLA_S2", "PLA_S3", "PLA_S4", "PLA_S5", 
                                      "BRE_S2", "BRE_S3", "BRE_S4", "BRE_S5","SAU_S2", "SAU_S3", "SAU_S4", "SAU_S5"))
## Graphe
ggplot(long_data, aes(x = Sample, y = Abundance, color = Sample)) +
  geom_jitter(width = 0.1) +
  facet_wrap(~ Pathway, scales = "free_y") +
  scale_color_manual(values = colors_sample) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing = unit(1, "lines"))
### Test connected VS disconnected
Carbon_cycle_abundance_coVSdeco <- as.data.frame(Carbon_cycle_abundance[!substr(rownames(Carbon_cycle_abundance), 1,3) %in% c("MOL", "PON"),]) # Removing MOL and PON, intermediate connectivity
Carbon_cycle_abundance_coVSdeco$Connectivity <- c(rep("Disconnected", 4), rep("Connected", 8), rep("Disconnected", 4),
                                                  rep("Disconnected", 4), rep("Disconnected", 4), rep("Connected", 4))
Carbon_cycle_abundance_coVSdeco$Sampling_campaign <- rep(c("May", "June", "August", "September"), 7)
Carbon_cycle_abundance_coVSdeco$Name <- substr(rownames(Carbon_cycle_abundance_coVSdeco), 1, 3)
## Generalized linear mixed models for each pathway
Carbon_cycle_abundance_coVSdeco_glm <- Carbon_cycle_abundance_coVSdeco[, -c(23:25)]
Carbon_pathways_coVSdeco <- colnames(Carbon_cycle_abundance_coVSdeco_glm)
Carbon_pathway_glmm <- lapply(Carbon_pathways_coVSdeco, function(pathway){
  # Temporary table
  data <- data.frame(
    Abondance = Carbon_cycle_abundance_coVSdeco_glm[[pathway]],
    Connectivity = Carbon_cycle_abundance_coVSdeco$Connectivity,
    Campaign = Carbon_cycle_abundance_coVSdeco$Sampling_campaign
  )
  # Model : 
  mod <- glmmTMB(
    Abondance ~ Connectivity + (1 | Campaign),
    family = nbinom2,
    data = data
  )
  # Null model
  mod_null <- glmmTMB(
    Abondance ~ (1 | Campaign),
    family = nbinom2,
    data = data
  )
  # Chi2 test
  lrt <- anova(mod_null, mod)
  p_global <- lrt$`Pr(>Chisq)`[2]
  chi2_global <- lrt$Chisq[2]
  # Post-hoc test
  posthoc <- emmeans(mod, pairwise ~ Connectivity)
  test <- as.data.frame(posthoc$contrasts)
  # Final table
  test$Pathway <- pathway
  test$Global_p <- p_global
  test$Chi2 <- chi2_global
  test
})
Carbon_pathway_stats <- do.call(rbind, Carbon_pathway_glmm)
unique(subset(Carbon_pathway_stats, Global_p < 0.05)$Pathway)
Carbon_cycle_abundance <- lapply(Carbon_cycle_abundance, as.numeric)
Carbon_cycle_abundance$Connectivity <- c(rep("Disconnected", 4), rep("Connected", 8), rep("Disconnected", 4), rep("Inter", 4), 
                                         rep("Disconnected", 4), rep("Inter", 4), rep("Disconnected", 4), rep("Connected", 4))
Carbon_cycle_abundance$Name <- substr(rownames(Carbon_cycle_abundance),1,3)
# Mean of pathways abundances co vs deco
Carbon_cycle_abundance_coVSdeco_moy <- Carbon_cycle_abundance_coVSdeco %>%
  group_by(Connectivity) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      )
    )
  )
Carbon_cycle_abundance_coVSdeco_moy <- Carbon_cycle_abundance_coVSdeco_moy[, -1]
rownames(Carbon_cycle_abundance_coVSdeco_moy) <- c("Lotic", "Lentic")

Carbon_cycle_abundance_len_sur_lot_ratio <- Carbon_cycle_abundance_coVSdeco_moy %>%
  summarise(
    across(
      ends_with("_mean"),
      ~ .x[rownames(Carbon_cycle_abundance_coVSdeco_moy) == "Lentic"] / .x[rownames(Carbon_cycle_abundance_coVSdeco_moy) == "Lotic"],
      .names = "{.col}_ratio"
    )
  )
hist(as.numeric(Carbon_cycle_abundance_len_sur_lot_ratio))

Carbon_cycle_abundance_lot_sur_len_ratio <- Carbon_cycle_abundance_coVSdeco_moy %>%
  summarise(
    across(
      ends_with("_mean"),
      ~ .x[rownames(Carbon_cycle_abundance_coVSdeco_moy) == "Lotic"] / .x[rownames(Carbon_cycle_abundance_coVSdeco_moy) == "Lentic"],
      .names = "{.col}_ratio"
    )
  )
hist(as.numeric(Carbon_cycle_abundance_lot_sur_len_ratio))
Carbon_cycle_abundance_ratio <- Carbon_cycle_abundance_len_sur_lot_ratio
Carbon_cycle_abundance_ratio[2,] <- Carbon_cycle_abundance_lot_sur_len_ratio
rownames(Carbon_cycle_abundance_ratio) <- c("Len/Lot", "Lot/Len")
  
Carbon_cycle_abundance_lone_coVSdeco_moy <- Carbon_cycle_abundance_coVSdeco %>%
  group_by(Name) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      )
    )
  )
Carbon_cycle_abundance_lone_coVSdeco_moy <- Carbon_cycle_abundance_lone_coVSdeco_moy[, -1]
rownames(Carbon_cycle_abundance_lone_coVSdeco_moy) <- unique(Carbon_cycle_abundance_coVSdeco$Name)
Carbon_cycle_abundance_coVSdeco$Name <- factor(Carbon_cycle_abundance_coVSdeco$Name, 
                                               levels = c("SAU", "BRE", "PLA", "MGA", "TON","MAT", "GAB"))



## Methane production and consumption
Methano_CO2_anova <- anova(lm(Methanogenesis_CO2 ~ Name, data = Carbon_cycle_abundance_coVSdeco, na.action = na.omit))
print(Methano_CO2_anova)
Methano_CO2_anova_aov <- aov(Methanogenesis_CO2 ~ Name, data = Carbon_cycle_abundance_coVSdeco)
Methano_CO2_anova_tukey <- TukeyHSD(Methano_CO2_anova_aov)
print(Methano_CO2_anova_tukey)

Methano_Acetate_anova <- anova(lm(Methanogenesis_Acetate ~ Name, data = Carbon_cycle_abundance_coVSdeco, na.action = na.omit))
print(Methano_Acetate_anova)
Methano_Acetate_anova_aov <- aov(Methanogenesis_Acetate ~ Name, data = Carbon_cycle_abundance_coVSdeco)
Methano_Acetate_anova_tukey <- TukeyHSD(Methano_Acetate_anova_aov)
print(Methano_Acetate_anova_tukey)

Methano_Methyle_anova <- anova(lm(Methanogenesis_Methyle ~ Name, data = Carbon_cycle_abundance_coVSdeco, na.action = na.omit))
print(Methano_Methyle_anova)
Methano_Methyle_anova_aov <- aov(Methanogenesis_Methyle ~ Name, data = Carbon_cycle_abundance_coVSdeco)
Methano_Methyle_anova_tukey <- TukeyHSD(Methano_Methyle_anova_aov)
print(Methano_Methyle_anova_tukey)

Methanotrophy_anova <- anova(lm(Methanotrophy_Methanol ~ Name, data = Carbon_cycle_abundance_coVSdeco, na.action = na.omit))
print(Methanotrophy_anova)
Methanotrophy_anova_aov <- aov(Methanotrophy_Methanol ~ Name, data = Carbon_cycle_abundance_coVSdeco)
Methanotrophy_anova_tukey <- TukeyHSD(Methanotrophy_anova_aov)
print(Methanotrophy_anova_tukey)

par(mfrow = c(2, 2))
boxplot(Methanogenesis_CO2 ~ Name, data = Carbon_cycle_abundance_coVSdeco, 
        col = c(rep("#D55E00", 4), rep("#009E73", 3)),
        ylim = c(0, 2000),
        xlab = "Sidearms",
        ylab = "Pathway Abundance (RPMs)", 
        main = "Hydrogenotrophic Methanogenesis")

boxplot(Methanogenesis_Acetate ~ Name, data = Carbon_cycle_abundance_coVSdeco, 
        col = c(rep("#D55E00", 4), rep("#009E73", 3)),
        ylim = c(0, 2000),
        xlab = "Sidearms",
        ylab = "Pathway Abundance (RPMs)", 
        main = "Aceticlastic Methanogenesis")

boxplot(Methanogenesis_Methyle ~ Name, data = Carbon_cycle_abundance_coVSdeco, 
        col = c(rep("#D55E00", 4), rep("#009E73", 3)),
        xlab = "Sidearms",
        ylab = "Pathway Abundance (RPMs)", 
        main = "Methylotrophic Methanogenesis")

boxplot(Methanotrophy_Methanol ~ Name, data = Carbon_cycle_abundance_coVSdeco, 
        col = c(rep("#D55E00", 4), rep("#009E73", 3)),
        xlab = "Sidearms",
        ylab = "Pathway Abundance (RPMs)", 
        main = "Methanotrophy")


###########