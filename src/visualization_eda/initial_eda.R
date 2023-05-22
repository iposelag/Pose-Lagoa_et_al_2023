#!/bin/Rscript
###############################################################################
#################### UPDATE ESET AND INITIAL EDA ##############################
###############################################################################
### This script has 2 parts:
## 1. Update the expression object with no negative expression genes
## 2. A first EDA analysis to understand data structure
## command: Rscript visualization_eda.R plt eset_name (plt: 14550 o 6480 & eset_name: expression)
## output: updated eset object and several pca, umap, clustering plots

## ----setup------------------------------------------------------------------------------------------------------------
# Set working directory
setwd("../../data/raw/GSE47460/")

# Set seet (reproducibility)
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(Biobase);library(cancerTiming); library(dplyr); library(factoextra)
library(dendextend); library(ggrepel); library(Rtsne); library(umap)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Set the platform of interest
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide two arguments.")
} else if (length(args)==2) {
  plt <- args[1]
  eset_name <- args[2]
}

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Create the required directories
if ("Visualization"%in%list.files("../../") == FALSE){dir.create("../../Visualization")}
if (plt%in%list.files("../../Visualization") == FALSE){
  dir.create(paste("../../Visualization/",plt, sep =""))}
# if ("all_samples"%in%list.files(paste("../../Visualization/", plt,"/",sep="")) == FALSE){
#   dir.create(paste("../../Visualization/",plt,"/all_samples", sep =""))}
# if ("filtered"%in%list.files(paste("../../Visualization/", plt,sep="")) == FALSE){
#   dir.create(paste("../../Visualization/",plt,"/filtered", sep =""))}

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Load data
load(paste("../../ExpressionSetObjects/",plt,"/",eset_name,".Rda", sep=""))
# Check data dimensions
dim(eset) # 17743x87
          # 18000x231 -> all_samples
          #
          # 18009x218 -> filtered
n_samples <- ncol(eset)
n_genes <- nrow(eset)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Deletion of negative expression genes
## ----------------------------------------------------------------------------------------------------------------------------------------
mypal <- colorRampPalette( c( "#3c556c", "#6c3c3d", "#556c3c" ) )(5)
negative_val <- rowSums(sign(exprs(eset))<0)>0
if(nrow(subset(exprs(eset),negative_val))>0){
  c <- rownames(subset(exprs(eset),negative_val))
  plot(exprs(eset)[c[1],],type= "l",ylim = c(0,25))
  for(i in 2:length(c)){
    lines(exprs(eset)[c[i],], col=mypal[i], lwd = 1)
  }
  eset <- eset[!(rownames(exprs(eset)) %in% c),]
  n_samples <- ncol(eset)
  n_genes <- nrow(eset)
}

# Check data dimensions
dim(eset)

# Update expression set (without genes negative expressions (if it corresponds))
save(eset,file = paste("../../ExpressionSetObjects/",plt,"/",eset_name,"_updated.Rda", sep=""))

## ----------------------------------------------------------------------------------------------------------------------------------------------
# EDA (inital to understand the data structure)
## ----------------------------------------------------------------------------------------------------------------------------------------------
# Outliers
if(plt == 6480){
  outliers <- c("GSM1149999", "GSM1150000", "GSM1150042", "GSM1150112") #6480 25%
}else{
  outliers <- c("GSM1150041", "GSM1150070", "GSM1150081", "GSM1150149", "GSM1150170", 
                "GSM1150187", "GSM1150246", "GSM1150413", "GSM1150499", "GSM1150502",
                "GSM1150517", "GSM1150523", "GSM1150537")}
## ----------------------------------------------------------------------------------------------------------------------------------------------
# Check gene expression levels
avgexp <- rowMeans(exprs(eset))
# png(paste("../../Visualization/",plt,"/all_samples/hist_exprs_levels.png", sep=""))
hist(avgexp, xlab = "expression levels", main = "Avge sample expression levels/gene", las=1)
# dev.off()

## ----------------------------------------------------------------------------------------------------------------------------------------------
# DataFrame of interest
exprs_values <- data.frame(t((Biobase::exprs(eset))))
exprs_values_ph <- exprs_values
exprs_values_ph$dis_condition <- as.factor(pData(eset)$dis_condition)
exprs_values_ph$sex <- as.factor(pData(eset)$sex)
exprs_values_ph$age <- cut(as.numeric(pData(eset)$age),
                           c((min(as.numeric(pData(eset)$age),na.rm=TRUE)-1),35,
                             45,55,65,75,max(as.numeric(pData(eset)$age),na.rm=TRUE)))
exprs_values_ph$smoker <- as.factor(pData(eset)$smoker)
exprs_values_ph$GOLD_stage <- as.factor(pData(eset)$GOLD_stage)
exprs_values_ph$pneumocystis_colonization <- as.factor(pData(eset)$pneumocystis_colonization)
exprs_values_ph$pred_dlco <- as.numeric(pData(eset)$pred_dlco)
exprs_values_ph$emphysema <- as.numeric(pData(eset)$emphysema)
exprs_values_ph$pred_fev_prebd <- as.numeric(pData(eset)$pred_fev_prebd)
exprs_values_ph$pred_fev_postbd <- as.numeric(pData(eset)$pred_fev_postbd)
exprs_values_ph$pred_fvc_prebd <- as.numeric(pData(eset)$pred_fvc_prebd)
exprs_values_ph$pred_fvc_postbd <- as.numeric(pData(eset)$pred_fvc_postbd)
exprs_values_ph$sample_id <- rownames(exprs_values_ph)
dim(exprs_values_ph)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Hierarchical Clustering

## Euclidean distance matrix
mat_dist <- get_dist(exprs_values, method="euclidean")

## Agglomeratie algorithm
# hc_complete <- hclust(d = mat_dist, method = "complete")
# hc_single <- hclust(d = mat_dist, method = "single")
hc_average <- hclust(d = mat_dist, method = "average")
hc_ward <- hclust(d = mat_dist, method = "ward.D2")
# hc_centroid <- hclust(d = mat_dist, method = "centroid")

p3 <- plot(x = hc_average, cex = 0.3, xlab = "", ylab = "", sub = "",
           main = "Linkage average")
p4 <- plot(x = hc_ward, cex = 0.3, xlab = "", ylab = "", sub = "",
           main = "Linkage ward")

## To know which linkage to use
cor(mat_dist, cophenetic(hc_average)) # 0.768
cor(mat_dist, cophenetic(hc_ward))

## WARD
hc_c_dendrogram <- as.dendrogram(hc_ward)

## Dendogram specifications
colors_to_use <- as.numeric(exprs_values_ph$dis_condition)
colors_to_use <- colors_to_use[order.dendrogram(hc_c_dendrogram)]
labels_colors(hc_c_dendrogram) <- colors_to_use
legend_colors <- seq(levels(exprs_values_ph$dis_condition))

# png(paste("../../Visualization/",plt,"/hc_averge.png",sep=""))
hc_c_dendrogram %>%
  set("labels_cex", 0.4) %>%
  plot(horiz= FALSE, axes=TRUE,
       main = "Cluster Dendrogram", ylab = "Euclidean Distance")
legend( "bottomright",
        cex=0.7,
        inset=c(0, 0.1),
        legend = levels(exprs_values_ph$dis_condition), fill=legend_colors)
# abline(h = 250 , col = 'purple')
rect.dendrogram(hc_c_dendrogram, k=2, border =1:2)
# dev.off()

## ----------------------------------------------------------------------------------------------------------------------------------------------
# PCA

pca <- prcomp(
  exprs_values,
  scale = TRUE,
  center = TRUE #we want the data scaled to have unit variance for each gene
)
head(pca$x[, 1:5])
pca_summary <- summary(pca)
pca_summary$importance[, 1:30]

## Make the first 30 PCs into a data frame for plotting with `ggplot2`
pca_df <- data.frame(pca$x[, 1:30]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::inner_join(
    dplyr::select(exprs_values_ph, sample_id, dis_condition, sex, age, smoker,
                  pneumocystis_colonization, GOLD_stage, emphysema, pred_fev_prebd,
                  pred_fev_postbd, pred_fvc_postbd, pred_fvc_prebd, pred_dlco),
    by = "sample_id"
  )


## Quality visualization
# png(paste("../../Visualization/",plt,"/quality_pca.png",sep=""))
fviz_eig(pca, n=15)
# dev.off()

## Eigenvalues
eig.val <- get_eigenvalue(pca)
eig.val[1:10,]

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
        scale_fill_continuous(type = "gradient",na.value="white")+
        geom_text_repel(data = subset(pca_df, sample_id %in% outliers),
                        aes(label = sample_id),
                        # family = "Poppins",
                        size = 2,
                        max.overlaps = Inf))
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
      color +
    geom_text_repel(data = subset(pca_df, sample_id %in% outliers),
                    aes(label = sample_id),
                    # family = "Poppins",
                    size = 2,
                    max.overlaps = Inf)
    )
  }
}


## Save the plots
# ggsave(paste("../../Visualization/",plt,"/all_samples/pca.png",sep=""), width = 15, height = 10)
# ggsave(paste("../../Visualization/",plt,"/filtered/pca.png",sep=""), width = 15, height = 10)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# K-means after PCA
## By Elbow method = total within-cluster sum of square (WSS)
fviz_nbclust(data.frame(pca$x[, 1:30]), kmeans, method = 'wss') +
  geom_vline(xintercept = 2, linetype = 2) +
  labs(subtitle = "Elbow method")
# ggsave(paste("../../Visualization/",plt,"/kmeans_elbow.png", sep=""), width = 15, height = 10)
## By Silhouette = for average silhouette width
fviz_nbclust(data.frame(pca$x[, 1:30]), kmeans, method = 'silhouette')
# ggsave(paste("../../Visualization/",plt,"/kmeans_silhouette.png", sep=""), width = 15, height = 10)
## By Gap Stats = for gap statistics
fviz_nbclust(data.frame(pca$x[, 1:30]), kmeans, method = 'gap_stat')
# ggsave(paste("../../Visualization/",plt,"/kmeans_gap_stat.png", sep=""), width = 15, height = 10)

k = 2
kmeans_pca = kmeans(data.frame(pca$x[, 1:30]), centers = k, nstart = 50)
fviz_cluster(kmeans_pca, data = data.frame(pca$x[, 1:30]))
## RESULTS
table(kmeans_pca$cluster, exprs_values_ph$dis_condition)


## ----------------------------------------------------------------------------------------------------------------------------------------------
# t-SNE
tSNE_fit <- exprs_values_ph[,1:n_genes] %>%
  scale() %>%
  Rtsne(perplexity=20)

## df
tSNE_df <- tSNE_fit$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  cbind(sample_id = exprs_values_ph[,"sample_id"]) %>%
  inner_join(exprs_values_ph[,(n_genes+1):ncol(exprs_values_ph)], by="sample_id")

## Plot
for(feature in c("emphysema", "pred_dlco", "pred_fev_postbd", 
                 "pred_fev_prebd", "pred_fvc_postbd", "pred_fvc_prebd", 
                 "dis_condition", "GOLD_stage", "sex", "smoker", "age", 
                 "pneumocystis_colonization")){ 
  if(is.numeric(pca_df[[feature]])){
    print(ggplot(
      tSNE_df,
      aes(
        x = tSNE1,
        y = tSNE2,
        shape = dis_condition,
        color = !!sym(feature) # label points with different colors for each `subgroup`
      )
    ) +
      geom_point() + # Plot individual points to make a scatterplot
      theme_classic() + # Format as a classic-looking plot with no gridlines
      scale_fill_continuous(type = "gradient",na.value="white")+
      geom_text_repel(data = subset(tSNE_df, sample_id %in% outliers),
                      aes(label = sample_id),
                      # family = "Poppins",
                      size = 2,
                      max.overlaps = Inf))
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
    tSNE_df,
    aes(
      x = tSNE1,
      y = tSNE2,
      shape = dis_condition,
      color = !!sym(feature) # label points with different colors for each `subgroup`
    )
  ) +
    geom_point() + # Plot individual points to make a scatterplot
    theme_classic() + # Format as a classic-looking plot with no gridlines
    color +
    geom_text_repel(data = subset(tSNE_df, sample_id %in% outliers),
                    aes(label = sample_id),
                    # family = "Poppins",
                    size = 2,
                    max.overlaps = Inf)
  )
}
## Save the plot
# ggsave(paste("../../Visualization/",plt,"/tSNE_plot.png", sep=""),height = 10, width = 15)


## ----------------------------------------------------------------------------------------------------------------------------------------------
# UMAP
umap_fit <- exprs_values_ph[1:n_genes] %>%
  scale() %>%
  umap()

## DF
umap_df <- umap_fit$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  cbind(sample_id = exprs_values_ph[,"sample_id"]) %>%
  inner_join(exprs_values_ph[,(n_genes+1):ncol(exprs_values_ph)], by="sample_id")

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
      scale_fill_continuous(type = "gradient",na.value="white")+
      geom_text_repel(data = subset(umap_df, sample_id %in% outliers),
                      aes(label = sample_id),
                      # family = "Poppins",
                      size = 2,
                      max.overlaps = Inf))
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
    color +
    geom_text_repel(data = subset(umap_df, sample_id %in% outliers),
                    aes(label = sample_id),
                    # family = "Poppins",
                    size = 2,
                    max.overlaps = Inf)
  )
}

## Save the plot
# ggsave(paste("../../Visualization/",plt,"/UMAP_plot.png",sep=""),height = 10, width = 15)


## ----------------------------------------------------------------------------------------------------------------------------------------------
# COPD VS GOLD STAGE
## ----------------------------------------------------------------------------------------------------------------------------------------

## DF
copd_gold <- filter(exprs_values_ph, dis_condition == "COPD")
copd_gold$GOLD_stage <- droplevels(copd_gold$GOLD_stage)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# PCA
pca <- prcomp(
  copd_gold[,1:n_genes],
  scale = TRUE,
  center = TRUE #we want the data scaled to have unit variance for each gene
)
head(pca$x[, 1:5])
pca_summary <- summary(pca)
pca_summary$importance[, 1:30]

## Make the first 30 PCs into a data frame for plotting with `ggplot2`
pca_df <- data.frame(pca$x[, 1:30]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::inner_join(
    dplyr::select(exprs_values_ph, sample_id, dis_condition, sex, age, smoker,
                  pneumocystis_colonization, GOLD_stage, emphysema, pred_fev_prebd,
                  pred_fev_postbd, pred_fvc_postbd, pred_fvc_prebd),
    by = "sample_id"
  )

## Quality visualization
# png(paste("../../Visualization/",plt,"/quality_pca_copd_gold.png",sep=""))
fviz_eig(pca)
# dev.off()

## Eigenvalues
# fviz_pca_var(pca, col.var="steelblue")
eig.val<-get_eigenvalue(pca)
eig.val[1:10,]

## Plot
ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = GOLD_stage # label points with different colors for each `subgroup`
  )
) +
  geom_point() + # Plot individual points to make a scatterplot
  theme_classic()  # Format as a classic-looking plot with no gridlines
  # geom_text_repel( data = subset(pca_df, sample_id %in% rownames(rare_samples)),
  #                  aes(label = sample_id),
  #                  family = "Poppins",
  #                  size = 2,
  #                  max.overlaps = Inf
  # )

## Save the plot
# ggsave(paste("../../Visualization/",plt,"/copd_gold.png",sep=""), width = 15, height = 10)
