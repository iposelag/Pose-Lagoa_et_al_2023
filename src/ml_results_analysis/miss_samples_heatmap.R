#!/bin/Rscript
###############################################################################
############## analysis of missclasified samples on cv and test ##############
###############################################################################


## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(sva); library(factoextra); library(dplyr); library(Biobase);
library(patchwork); library(pheatmap); library(purrr); library(viridis)
library(tidyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# LOAD eset object
load("../../data/ExpressionSetObjects/eset_deg_lfc.Rda")
dim(eset_deg)
n_samples <- ncol(eset_deg)
n_genes <- nrow(eset_deg)

methodology <- "optm/disgenet_expansion/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories


for (name in c("train", "test")) {

  # Load miss_samples table
  miss_samples <- read.csv(paste("../../data/ML_out/two_class/",methodology,"miss_samples/miss_samples_",name,"_pheno.csv", sep=""))
  rownames(miss_samples) <- miss_samples$sample_id
  # Save intersection samples for further analysis
  assign(paste("miss_",name,sep=""),miss_samples)
  assign(paste("miss_intersection_",name,sep=""),miss_samples[complete.cases(miss_samples[, 2:7]), ])
  # miss_samples <- get(paste("miss_intersection_",name,sep=""))
  
  pheno <- miss_samples[,c(8:length(miss_samples))]
  # pheno$sample_id <- miss_samples$sample_id
  miss_samples_plot <- miss_samples[,c(2:7)]
  miss_samples_plot[is.na(miss_samples_plot)] <- 0

  #Save order or pheatmap clustering
  pheatm <- pheatmap(as.matrix(miss_samples_plot), cluster_cols = FALSE)
  miss_samples_plot$sample_id <- miss_samples$sample_id
  miss_samples_plot <- cbind(miss_samples_plot, pheno)
  # miss_samples_plot <- merge(miss_samples_plot, pheno, by="sample_id")
  samples_order <- miss_samples_plot[pheatm$tree_row$order,]$sample_id
  
  # Reshape data from wide to long format
  df_long <- gather(miss_samples_plot, key = "metric", value = "value", -sample_id, -dis_condition, -GOLD_stage, -sex, -smoker,
                    - age, -pneumocystis_colonization, -emphysema, -pred_dlco, -pred_fev_prebd, -pred_fev_postbd,
                    -pred_fvc_postbd, -pred_fvc_prebd, -platform_id)
  df_long$value <- as.numeric(df_long$value)

  # Order sample_id based on dis_condition (needed to have our heatmap ordered by condition)
  # df_long <- df_long %>% arrange(dis_condition)
  # df_long$sample_id <- factor(df_long$sample_id, levels = unique(df_long[order(df_long$dis_condition), "sample_id"]))

  # Create heatmap
  heatmap <- ggplot(df_long, aes(x = sample_id, y = metric, fill = value)) +
    geom_tile(color = "white",
               linetype = 1,
               width = 0.8,
               height = 0.8) +
    scale_fill_gradient(low = "#C1D9D0", high = "#669BF0", na.value = "white")+
    theme(text = element_text(size=8),
          axis.text.x = element_text(angle=90, hjust=1))+
    xlim(samples_order)


  for(feature in c("emphysema", "pred_dlco", "pred_fev_postbd",
                   "pred_fev_prebd", "pred_fvc_postbd", "pred_fvc_prebd",
                   "dis_condition", "GOLD_stage", "sex", "smoker", "age",
                   "pneumocystis_colonization")){
    
    if(is.numeric(df_long[[feature]])){
      assign(paste("side_bar_",feature, sep =""), ggplot(df_long, aes(x = sample_id, y = "", fill = !!sym(feature))) +
               geom_col(width = 1) +
               theme_void() +
               scale_fill_continuous(type = "viridis",na.value="white")+
               theme(legend.position = "right",
                     legend.direction = "horizontal",
                     legend.key.size = unit(0.3, 'cm'),
                     legend.text = element_text(size=6),
                     legend.title = element_text(size=8),
                     axis.title = element_blank())+ 
               xlim(samples_order) )
      next
    }else if (feature == "dis_condition"){
      color <- scale_fill_manual(values = c("CTRL" = "#ffdcdb", "COPD" = "#91a8d0"),
                                 na.value = "white")
    }else if (feature == "GOLD_stage") {
      color <- scale_fill_manual(values = c("0-At Risk" = "#f6e0b5", "1-Mild COPD" = "#eea990",
                                            "2-Moderate COPD" = "#aa6f73", "3-Severe COPD" = "#a39193",
                                            "4-Very Severe COPD" = "#66545e"),
                                 na.value = "white")
    }else if (feature == "sex") {
      color <- scale_fill_manual(values = c("1-Male" = "#e1f7d5", "2-Female" = "#c9c9ff"),
                                 na.value = "white")
    }else if (feature== "smoker") {
      color <- scale_fill_manual(values = c("1-Current" = "#83adb5", "2-Ever (>100)" = "#c7bbc9",
                                            "3-Never"="#5e3c58"),
                                 na.value = "white")
    }else if (feature == "age"){
      color <- scale_fill_manual(values = c("(27,35]" = "#ece6ff", "(35,45]" = "#efbbff",
                                            "(45,55]" = "#d896ff", "(55,65]" = "#be29ec",
                                            "(65,75]" = "#800080", "(75,91]" = "#660066"),
                                 na.value = "white")

    }else if (feature == "pneumocystis_colonization"){
      color <- scale_fill_manual(values = c("Negative" = "#ffbaba",
                                            "Positive" ="#b8d8be"),
                                 na.value = "white")

    }
    assign(paste("side_bar_",feature, sep =""), ggplot(df_long, aes(x = sample_id, y = "", fill = !!sym(feature))) +
             geom_col(width = 1) +
             theme_void() +
             color +
             theme(legend.position = "right",
                   legend.direction = "horizontal",
                   legend.key.size = unit(0.3, 'cm'),
                   legend.text = element_text(size=6),
                   legend.title = element_text(size=8),
                   axis.title = element_blank())+
                   #  axis.text.x = element_text(angle = 90, hjust = 1))+ 
             xlim(samples_order))


  }
  plot <- heatmap +
    side_bar_dis_condition + side_bar_GOLD_stage + side_bar_sex + side_bar_smoker +
    side_bar_pneumocystis_colonization +  side_bar_age +
    side_bar_emphysema + side_bar_pred_fev_prebd +side_bar_pred_fev_postbd +
    side_bar_pred_fvc_prebd + side_bar_pred_fvc_postbd + side_bar_pred_dlco +
    plot_layout(ncol = 1, heights = c(2,0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1))
  # plot <- heatmap + side_bar_dis_condition + side_bar_GOLD_stage + side_bar_sex + side_bar_age +
  #   side_bar_emphysema + side_bar_pred_fev_prebd +side_bar_pred_fev_postbd +
  #  side_bar_pred_dlco +side_bar_pred_fvc_prebd +side_bar_pred_fvc_postbd+
  #   plot_layout(ncol = 1, heights = c(2,0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1))

  ggsave(paste("../../data/ML_out/two_class/",methodology,"miss_samples/histogram_",name,".pdf",sep=""),
         width = 12, height = 10)


}



# Enrichment analysis

pdata <- miss_train
head(pdata)
pdata$miss_sample <- ifelse(is.na(miss_train$Freq_svm_p), 0, 1)
pdata <- pdata[,-c(2:7)]
head(pdata)

pdata <- pdata %>%
  mutate(dis_condition = factor(dis_condition),
         sex = factor(sex),
         pneumocystis_colonization = factor(pneumocystis_colonization),
         smoker = factor(smoker),
         GOLD_stage = factor(GOLD_stage),
         platform_id = factor(platform_id),
         age = as.numeric(age),
         pred_dlco = as.numeric(pred_dlco),
         pred_fev_prebd = as.numeric(pred_fev_prebd),
         pred_fev_postbd = as.numeric(pred_fev_postbd),
         pred_fvc_prebd = as.numeric(pred_fvc_prebd),
         pred_fvc_postbd = as.numeric(pred_fvc_postbd))

for(i in 2:(ncol(pdata)-1)){
  print(names(pdata[i]))
    if(class(pdata[[i]]) != "factor"){
      print(ggplot(pdata, aes(x=pdata[[i]], fill=as.factor(miss_sample))) +
        geom_density(alpha=0.5) +
        xlab(names(pdata[i])) +
        ylab("Density") +
        scale_fill_manual(values=c("blue", "green"), name="Missclasified Samples"))
      normality <- by(pdata[[i]], pdata$miss_sample, shapiro.test)
      if(normality$`0`$p.value < 0.05 | normality$`1`$p.value < 0.05){
        print("Variable does not follow normality:")
        if(wilcox.test(pdata[[i]]~pdata$miss_sample)$p.value < 0.05){
          print("The variable is enriched in missclasified samples:")
          print(wilcox.test(pdata[[i]]~pdata$miss_sample))
        }else{
          print("The variable is not enriched in missclasified samples")
        }
      }else if(var.test(pdata[[i]]~pdata$miss_sample)$p.value < 0.05){
        print("Variable does not follow homocedasticity:")
        if(wilcox.test(pdata[[i]]~pdata$miss_sample)$p.value < 0.05){
          print("The variable is enriched in missclasified samples:")
          print(wilcox.test(pdata[[i]]~pdata$miss_sample))
        }else{
          print("The variable is not enriched in missclasified samples")
        }

      }else{
        print("Variable follows normality and homocedasticity")
        if(t.test(pdata[[i]]~pdata$miss_sample)$p.value < 0.05){
          print("The variable is enriched in missclasified samples:")
          print(t.test(pdata[[i]]~pdata$miss_sample))
        }else{
          print("The variable is not enriched in missclasified samples")
        }

      }

    }else{
      tab <- table(pdata$miss_sample, pdata[[i]])
      # ggplot(pdata, aes(x = miss_sample, y = pdata[[i]], fill = rownames(pdata[[i]]))) +
      #   geom_bar() +
      #   scale_fill_manual(values = c("#4e79a7", "#f28e2c"),
      #                     name = "Row",
      #                     labels = c("0", "1")) +
      #   labs(x = "Platform", y = "Count", title = "Bar plot of GPL14550 and GPL6480") +
      #   theme_minimal()
      if(length(levels(pdata[[i]])) == 2){
        if(chisq.test(tab, correct = TRUE)$p.value < 0.05){
          print("The variable is enriched in missclasified samples:")
          test <- chisq.test(tab, correct = TRUE)
          print(test)
          print(test$expected); print(test$observed)
        }else{
          print("The variable is not enriched in missclasified samples")
        }
      }else {
        if(chisq.test(tab, simulate.p.value = TRUE)$p.value < 0.05){
          print("The variable is enriched in missclasified samples:")
          test <- chisq.test(tab, simulate.p.value = TRUE)
          print(test)
          print(test$expected); print(test$observed)
        }else{
          print("The variable is not enriched in missclasified samples")
        }
    }
  }
}


ggplot(proba) +
  aes(
    x = pdata.GOLD_stage,
    fill = pdata.GOLD_stage,
    weight = as.factor(pdata.miss_sample)
  ) +
  geom_bar() +
  scale_fill_hue() +
  theme_minimal() +
  facet_wrap(vars(as.factor(pdata$miss_sample)))

ggplot(data_to_plot) +
  aes(
    x = methodology,
    fill = classifier,
    group = classifier,
    weight = estimate
  ) +
  labs(x = "accuracy")+
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(rf = "#6b3e26",
               svm_r = "#ffc5d9",
               svm_p = "#c2f2d0",
               knn = "#ffcb85",
               pen_reg = "#fdf5c9",
               xgb = "#ff6f69"
    )
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Performance bewteen classifiers and methodologies",
    subtitle = "Train"
  ) +
  geom_text(
    aes(label = round(estimate, 2), y =estimate),
    hjust = -0.5,
    size = 2,
    position = position_dodge(width = 0.9)
  ) +
  coord_flip() +
  theme_minimal()

# NA
# mellorar plot
# plot tab

pdata$emphysema <- replace(pdata$emphysema, is.na(pdata$emphysema), 0)

ggplot(pdata, aes(x=emphysema, fill=as.factor(miss_sample))) +
  geom_density(alpha=0.5) +
  ggtitle("Distribution of emphysema for both groups") +
  xlab("Emphysema") +
  ylab("Density") +
  scale_fill_manual(values=c("blue", "green"), name="Missclasified Samples")

# write.csv(miss_intersection_methods,paste( "../../data/ML_out/two_class/",methodology,"/dea/miss_samples/miss_intersection_methods.csv",sep=""))

# PCA
pca_copd_ctrl <- prcomp(
  combat_edata_ph[,1:nrow(combat_edata)],
  scale = TRUE, # we want the data scaled to have unit variance for each gene
  center = TRUE
)
head(pca_copd_ctrl$x[, 1:5])
pca_copd_ctrl_summary <- summary(pca_copd_ctrl)
pca_copd_ctrl_summary$importance[, 1:10]
fviz_eig(pca_copd_ctrl)

pca_df_copd_ctrl <- data.frame(pca_copd_ctrl$x[, 1:10]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::inner_join(
    dplyr::select(combat_edata_ph, sample_id, dis_condition, age, sex, GOLD_stage,
                  smoker, pneumocystis_colonization, emphysema, pred_dlco, platform_id),
    by = "sample_id"
  )


## Eigenvalues
eig.val<-get_eigenvalue(pca_copd_ctrl)
eig.val[1:10,]

library(ggrepel)
library(viridis)
p2 <- ggplot(
  pca_df_copd_ctrl,
  aes(
    x = PC1,
    y = PC2,
    shape = dis_condition,
    color = platform_id # label points with different colors for each `subgroup`
  ))  +
  geom_point() + # Plot individual points to make a scatterplot
  labs(x = "PC1",
       y = "PC2",
       subtitle = "PCA plot of deg from COPD+CTRL comparison (combat matrix) ") +
  theme_classic() + # Format as a classic-looking plot with no gridlines
  #scale_color_manual(values=c('black','#f2a297','#b0d1b2')) +
  #scale_color_viridis() +
  geom_text_repel( data = filter(pca_df_copd_ctrl, sample_id %in%
                                   miss_intersection_methods$sample_id),
                   aes(label = sample_id),
                   family = "Poppins",
                   size = 2,
                   max.overlaps = Inf
  )

p2
ggsave(paste("../../data/ML_out/two_class/",methodology,"/dea/miss_samples/pca_miss_intersection.png",sep=""), width = 12, height = 8)

# genes set analysis
# esquema prefiltering con catidades q vou perdendo
# plots

# pData(eset) <- rename(pData(eset), sample_id = geo_accession)
# head(miss_samples)
# miss_samples_res[miss_samples_res == TRUE] <- 'train'
# miss_samples_res <- rename(miss_samples_res,
#                          RF = Freq_RF,
#                          SVM_R = Freq_SVM_r,
#                          SVM_P = Freq_SVM_p,
#                          PEN_REG = Freq_pen_reg,
#                          KNN = Freq_knn,
#                          XGB = Freq_xgb)
# miss_all_train[miss_all_train == TRUE] <- "train_all"
# miss_test[miss_test == TRUE] <- "test"
#
#
# pData(eset) <- merge(pData(eset), miss_samples_res, by = "sample_id", all =TRUE)
# pData(eset)  <- merge(pData(eset), miss_all_train,  by = "sample_id", all =TRUE)
# pData(eset)  <- merge(pData(eset), miss_test, by = "sample_id", all = TRUE)
#
# save(eset, file = "../../data/ExpressionSetObjects/eset_jon.Rda")
