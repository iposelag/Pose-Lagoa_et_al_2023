#!/bin/Rscript
###############################################################################
################# OBTAIN ESET OBJECTS BOTH PLATFORM ###########################
###############################################################################
## output: new eset object (with expression data from both platforms) and several pca and umap plots

# Set seed reproducibility
set.seed(1234)
## This option avoid use_directory() being verbose later on

## --------------------------------------------------------------------------------------------------------------------
# Required packages
library(cancerTiming); library(limma); library(sva); library(ggplot2); library(ggrepel)
library(tidyr); library(factoextra); library(Biobase); library(msigdbr)
library(clusterProfiler); library(dendextend); library("ReactomePA")
require(DOSE);library(enrichplot); library(umap); library(dplyr)

## --------------------------------------------------------------------------------------------------------------------
# Create required directories
# if ("RESULTS"%in%list.files(".") == FALSE){dir.create("RESULTS", showWarnings = FALSE)}
# if ("H_CLUST"%in%list.files("RESULTS") == FALSE){dir.create("RESULTS/H_CLUST", showWarnings = FALSE)}

if ("both_plt"%in%list.files("../../Visualization/") == FALSE){dir.create("../../Visualization/both_plt", showWarnings = FALSE)}

## --------------------------------------------------------------------------------------------------------------------
# Set working directory
setwd("../../data/Visualization/both_plt/")

## --------------------------------------------------------------------------------------------------------------------
# Upload data

# Platform 6480
load("../../ExpressionSetObjects/6480/expression_updated.Rda")
eset_6480 <- eset
print("Dimensions eset_6480")
dim(eset_6480)  # 17743x87 all samples
# 17781x83 filtered
# Platform 14550
load("../../ExpressionSetObjects/14550/expression_updated.Rda")
eset_14550 <- eset
print("Dimensions eset_14550")
dim(eset_14550)  # 17995x231 all samples
# 18004x218 filtered

## --------------------------------------------------------------------------------------------------------------------
# Combine platforms
## genes intersection
genes_intersection <- Reduce(intersect, list(fData(eset_14550)$GENE_SYMBOL, fData(eset_6480)$GENE_SYMBOL))

## plt 14550
fData(eset_14550) <- filter(fData(eset_14550), GENE_SYMBOL %in% genes_intersection)
eset_14550 <- eset_14550[fData(eset_14550)$GENE_SYMBOL,]
table(rownames(eset_14550) == fData(eset_14550)$GENE_SYMBOL)
print("Dimensions eset_14550: intersection genes")
dim(eset_14550)
samples_14550 <- eset_14550@phenoData@data[["geo_accession"]]
## Save the data
save(samples_14550,file = "../../ExpressionSetObjects/14550/samples_14550.Rda")

## plt 6480
fData(eset_6480) <- filter(fData(eset_6480), GENE_SYMBOL %in% genes_intersection)
eset_6480 <- eset_6480[fData(eset_6480)$GENE_SYMBOL,]
table(rownames(eset_6480) == fData(eset_6480)$GENE_SYMBOL)
print("Dimensions eset_6480: intersection genes")
dim(eset_6480)
samples_6480 <- eset_6480@phenoData@data[["geo_accession"]]
## Save the data
save(samples_6480,file = "../../ExpressionSetObjects/6480/samples_6480.Rda")


# Combination
eset <- BiocGenerics::combine(eset_14550,eset_6480)
# Check data dimensions
print("Dimensions eset: platforms combination")
dim(eset) # 16230x318 all_Samples
# 16235x301 filtered
n_samples <- ncol(eset)
n_genes <- nrow(eset)

# Save the data
save(eset, file = "../../ExpressionSetObjects/eset.Rda")
# load("../../data/ExpressionSetObjects/eset_deg_lfc.Rda")
# eset <- eset_deg

## --------------------------------------------------------------------------------------------------------------------
# EDA both platforms
## --------------------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------------------
# BATCH effect (platform)

## --------------------------------------------------------------------------------------------------------------------
# Batch effect detection (also see the results of AQM)
# png("Plots/mds.png", heigh=1000, width = 1000 )
plotMDS(eset)
# dev.off()
# a<-"lfc"
# load(paste("../../ExpressionSetObjects/eset_deg_",a,".Rda",sep=""))
# eset <- eset_deg
## --------------------------------------------------------------------------------------------------------------------
# BATCH effect correction (platform): ComBAT
batch = pData(eset)$platform_id
#' model w effect of interest
modcombat = model.matrix(~ dis_condition, data=pData(eset))
combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)


## --------------------------------------------------------------------------------------------------------------------
# Data frame of interest

combat_edata_ph <- as.data.frame(t(combat_edata))
combat_edata_ph$dis_condition <- as.factor(pData(eset)$dis_condition)
combat_edata_ph$sex <- as.factor(pData(eset)$sex)
combat_edata_ph$age <- cut(as.numeric(pData(eset)$age),
                           c((min(as.numeric(pData(eset)$age),na.rm=TRUE)-1),35,
                             45,55,65,75,max(as.numeric(pData(eset)$age),na.rm=TRUE)))
combat_edata_ph$smoker <- as.factor(pData(eset)$smoker)
combat_edata_ph$GOLD_stage <- as.factor(pData(eset)$GOLD_stage)
combat_edata_ph$pneumocystis_colonization <- as.factor(pData(eset)$pneumocystis_colonization)
combat_edata_ph$pred_dlco <- as.numeric(pData(eset)$pred_dlco)
combat_edata_ph$emphysema <- as.numeric(pData(eset)$emphysema)
combat_edata_ph$pred_fev_prebd <- as.numeric(pData(eset)$pred_fev_prebd)
combat_edata_ph$pred_fev_postbd <- as.numeric(pData(eset)$pred_fev_postbd)
combat_edata_ph$pred_fvc_prebd <- as.numeric(pData(eset)$pred_fvc_prebd)
combat_edata_ph$pred_fvc_postbd <- as.numeric(pData(eset)$pred_fvc_postbd)
combat_edata_ph$sample_id <- rownames(combat_edata_ph)

## --------------------------------------------------------------------------------------------------------------------
# PCA COPD+CTRL
print("PCA analysis")
pca <- prcomp(
  data.frame(t(combat_edata)),
  scale = TRUE,
  center = TRUE #we want the data scaled to have unit variance for each gene
)
head(pca$x[, 1:5])
pca_summary <- summary(pca)

# Importance of the first 30 pc
pca_summary$importance[, 1:30]

fviz_eig(pca, ncp = 15)

# First 30 PCs into a data frame
pca_df <- data.frame(pca$x[, 1:30]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::inner_join(
    dplyr::select(combat_edata_ph, sample_id, dis_condition, sex, age, smoker,
                  pneumocystis_colonization, GOLD_stage, pred_dlco, emphysema, pred_fev_prebd, 
                  pred_fev_postbd, pred_fvc_prebd, pred_fvc_postbd),
    by = "sample_id"
  )

# PCA plot
## Plot
for(pc in c("PC1","PC3")){
  for(feature in c("emphysema", "pred_dlco", "pred_fev_postbd", 
                   "pred_fev_prebd", "pred_fvc_postbd", "pred_fvc_prebd", 
                   "dis_condition", "GOLD_stage", "sex", "smoker", "age", 
                   "pneumocystis_colonization")){ 
    if(is.numeric(pca_df[[feature]])){
      print(ggplot(
        pca_df,
        aes(
          x = !!sym(pc),
          y = PC2,
          shape = dis_condition,
          color = !!sym(feature) # label points with different colors for each `subgroup`
        )
      ) +
        geom_point() + # Plot individual points to make a scatterplot
        theme_classic() + # Format as a classic-looking plot with no gridlines
        scale_fill_continuous(type = "gradient",na.value="white"))
      # geom_text_repel(data = subset(pca_df, sample_id %in% outliers),
      #                 aes(label = sample_id),
      #                 # family = "Poppins",
      #                 size = 2,
      #                 max.overlaps = Inf))
      next
    }else if (feature == "dis_condition"){
      color <- scale_color_manual(values = c("CTRL" = "#f7cac9", "COPD" = "#91a8d0"))
    }else if (feature == "GOLD_stage") {
      color <- scale_color_manual(values = c("0-At Risk" = "#f6e0b5", "1-Mild COPD" = "#eea990", 
                                             "2-Moderate COPD" = "#aa6f73", "3-Severe COPD" = "#a39193",
                                             "4-Very Severe COPD" = "#66545e"))
    }else if (feature == "sex") {
      color <- scale_color_manual(values = c("1-Male" = "#b8d8be", "2-Female" = "#c9c9ff"))
    }else if (feature== "smoker") {
      color <- scale_color_manual(values = c("1-Current" = "#83adb5", "2-Ever (>100)" = "#c7bbc9",
                                             "3-Never"="#5e3c58"))
    }else if (feature == "age"){
      color <- scale_color_manual(values = c("(27,35]" = "#ece6ff", "(35,45]" = "#efbbff", 
                                             "(45,55]" = "#d896ff", "(55,65]" = "#be29ec", 
                                             "(65,75]" = "#800080", "(75,91]" = "#660066"))
      
    }else if (feature == "pneumocystis_colonization"){
      color <- scale_color_manual(values = c("Negative" = "#ffbaba",
                                             "Positive" ="#b8d8be"))
      
    }
    print(ggplot(
      pca_df,
      aes(
        x = !!sym(pc),
        y = PC2,
        shape = dis_condition,
        color = !!sym(feature) # label points with different colors for each `subgroup`
      )
    ) +
      geom_point() + # Plot individual points to make a scatterplot
      theme_classic() + # Format as a classic-looking plot with no gridlines
      color 
    # geom_text_repel(data = subset(pca_df, sample_id %in% outliers),
    #                 aes(label = sample_id),
    #                 # family = "Poppins",
    #                 size = 2,
    #                 max.overlaps = Inf)
    )
  }
}

## ----------------------------------------------------------------------------------------------------------------------------------------------
# UMAP
umap_fit <- combat_edata_ph[1:n_genes] %>%
  scale() %>%
  umap()

## DF
umap_df <- umap_fit$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  cbind(sample_id = combat_edata_ph[,"sample_id"]) %>%
  inner_join(combat_edata_ph[,(n_genes+1):ncol(combat_edata_ph)], by="sample_id")

## PLOT
for(feature in c("emphysema", "pred_dlco", "pred_fev_postbd", 
                 "pred_fev_prebd", "pred_fvc_postbd", "pred_fvc_prebd", 
                 "dis_condition", "GOLD_stage", "sex", "smoker", "age", 
                 "pneumocystis_colonization")){ 
  if(is.numeric(pca_df[[feature]])){
    print(ggplot(
      umap_df,
      aes(
        x = UMAP1,
        y = UMAP2,
        shape = dis_condition,
        color = !!sym(feature) # label points with different colors for each `subgroup`
      )
    ) +
      geom_point() + # Plot individual points to make a scatterplot
      theme_classic() + # Format as a classic-looking plot with no gridlines
      scale_fill_continuous(type = "gradient",na.value="white"))
    # geom_text_repel(data = subset(umap_df, sample_id %in% outliers),
    #                 aes(label = sample_id),
    #                 # family = "Poppins",
    #                 size = 2,
    #                 max.overlaps = Inf))
    next
  }else if (feature == "dis_condition"){
    color <- scale_color_manual(values = c("CTRL" = "#f7cac9", "COPD" = "#91a8d0"))
  }else if (feature == "GOLD_stage") {
    color <- scale_color_manual(values = c("0-At Risk" = "#f6e0b5", "1-Mild COPD" = "#eea990", 
                                           "2-Moderate COPD" = "#aa6f73", "3-Severe COPD" = "#a39193",
                                           "4-Very Severe COPD" = "#66545e"))
  }else if (feature == "sex") {
    color <- scale_color_manual(values = c("1-Male" = "#b8d8be", "2-Female" = "#c9c9ff"))
  }else if (feature== "smoker") {
    color <- scale_color_manual(values = c("1-Current" = "#83adb5", "2-Ever (>100)" = "#c7bbc9",
                                           "3-Never"="#5e3c58"))
  }else if (feature == "age"){
    color <- scale_color_manual(values = c("(27,35]" = "#ece6ff", "(35,45]" = "#efbbff", 
                                           "(45,55]" = "#d896ff", "(55,65]" = "#be29ec", 
                                           "(65,75]" = "#800080", "(75,91]" = "#660066"))
    
  }else if (feature == "pneumocystis_colonization"){
    color <- scale_color_manual(values = c("Negative" = "#ffbaba",
                                           "Positive" ="#b8d8be"))
    
  }
  print(ggplot(
    umap_df,
    aes(
      x = UMAP1,
      y = UMAP2,
      shape = dis_condition,
      color = !!sym(feature) # label points with different colors for each `subgroup`
    )
  ) +
    geom_point() + # Plot individual points to make a scatterplot
    theme_classic() + # Format as a classic-looking plot with no gridlines
    color 
  # geom_text_repel(data = subset(umap_df, sample_id %in% outliers),
  #                 aes(label = sample_id),
  #                 # family = "Poppins",
  #                 size = 2,
  #                 max.overlaps = Inf)
  )
}

## --------------------------------------------------------------------------------------------------------------------
# Generate an ESET with 3 class: CTRL, SEV, MSEV
# load(file = "../../data/ExpressionSetObjects/eset.Rda")
# eset_several_class <- eset

# eset_several_class$GOLD_stage <- as.factor(pData(eset_several_class)$GOLD_stage)
# eset_several_class$GOLD_stage <- recode(eset_several_class$GOLD_stage, "0-At Risk"="At_Risk",
#                                      "1-Mild COPD"="No_Severe", "2-Moderate COPD"="No_Severe",
#                                      "3-Severe COPD"="Severe", "4-Very Severe COPD"="Severe")

# eset_several_class<- eset_several_class[,!is.na(eset_several_class$GOLD_stage)]
# dim(eset_several_class)
# eset_several_class$dis_condition <- as.factor(paste(eset_several_class$dis_condition, eset_several_class$GOLD_stage, sep="_"))
# eset_several_class$dis_condition <- relevel(eset_several_class$dis_condition,
#                                           ref="COPD_Severe")
# save(eset_several_class, file = "../../data/ExpressionSetObjects/eset_several_class.Rda")
