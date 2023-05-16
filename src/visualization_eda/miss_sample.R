## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(sva); library(factoextra); library(dplyr); library(Biobase); 
library(UpSetR); library(patchwork); library(pheatmap); library(purrr);
library(tidyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# LOAD eset object
load("../../data/ExpressionSetObjects/eset_deg_lfc.Rda")
dim(eset_deg)
n_samples <- ncol(eset_deg)
n_genes <- nrow(eset_deg)

methodology <- "optm/data_driven/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("miss_samples"%in%list.files(paste("../../data/ML_out_out/two_class/",methodology, sep="")) == FALSE){
  dir.create(paste("../../data/ML_out/two_class/",methodology,"miss_samples", sep=""), showWarnings = FALSE)}


# Load ML results
load(paste("../../data/ML_out/two_class/",methodology, "results.Rda", sep= ""))

batch = pData(eset_deg)$platform_id
modcombat = model.matrix(~ dis_condition, data=pData(eset_deg))
combat_edata = ComBat(dat=exprs(eset_deg), batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)
combat_edata_ph <- as.data.frame(t(combat_edata))
combat_edata_ph$dis_condition <- as.factor(pData(eset_deg)$dis_condition)
combat_edata_ph$sex <- as.factor(pData(eset_deg)$sex)
combat_edata_ph$age <- as.numeric(pData(eset_deg)$age)
# combat_edata_ph$age <- cut(as.numeric(pData(eset_deg)$age),
#                            c((min(as.numeric(pData(eset_deg)$age),na.rm=TRUE)-1),35,45,55,65,75,max(as.numeric(pData(eset_deg)$age),na.rm=TRUE)))
combat_edata_ph$smoker <- as.factor(pData(eset_deg)$smoker)
combat_edata_ph$GOLD_stage <- as.factor(pData(eset_deg)$GOLD_stage)
combat_edata_ph$pneumocystis_colonization <- as.factor(pData(eset_deg)$pneumocystis_colonization)
combat_edata_ph$emphysema <- as.numeric(pData(eset_deg)$emphysema)
combat_edata_ph$pred_dlco <- as.numeric(pData(eset_deg)$pred_dlco)
# combat_edata_ph$pred_dlco <- cut(as.numeric(pData(eset_deg)$pred_dlco),
#                                  c((min(as.numeric(pData(eset_deg)$pred_dlco),na.rm=TRUE)-1),39,60,75,max(as.numeric(pData(eset_deg)$pred_dlco),na.rm=TRUE)))
# levels(combat_edata_ph$pred_dlco) <- c("severe", "moderate", "mild", "normal")
combat_edata_ph$pred_fev_prebd <- as.numeric(pData(eset_deg)$pred_fev_prebd)
combat_edata_ph$pred_fev_postbd <- as.numeric(pData(eset_deg)$pred_fev_postbd)
combat_edata_ph$pred_fvc_prebd <- as.numeric(pData(eset_deg)$pred_fvc_prebd)
combat_edata_ph$pred_fvc_postbd <- as.numeric(pData(eset_deg)$pred_fvc_postbd)
combat_edata_ph$dis_condition <- relevel(combat_edata_ph$dis_condition,
                                         ref="COPD")
combat_edata_ph$sample_id <- rownames(combat_edata_ph)
combat_edata_ph$platform_id<- as.factor(pData(eset_deg)$platform_id)
dim(combat_edata_ph)

pheno_cctics <- c("dis_condition", "sample_id","GOLD_stage",
                  "sex", "smoker", "age", "pneumocystis_colonization",
                  "emphysema", "pred_dlco", "pred_fev_prebd", "pred_fev_postbd",
                  "pred_fvc_prebd", "pred_fvc_postbd", "platform_id")

splits <- c("miss_samples_train_res", "miss_samples_all_train", "miss_samples_test")
splt <- c("train_res", "all_train", "test")
classifiers <- c("rf","svm_r","svm_p","pen_reg","knn", "xgb")

for(j in 1:length(splits)){
  split <- splits[j]
  for(i in 1:length(classifiers)){
    
    classifier <- classifiers[i]
    assign(paste("miss_",classifier,"_",splt[j], sep =""),results[[splits[j]]][[i]])
    
  }
}

# Missclasified train res add samples as a column and delete row info
## rf
miss_rf_train_res$sample_id <- rownames(miss_rf_train_res); miss_rf_train_res <- miss_rf_train_res[,-c(1)]
## svm_r
miss_svm_r_train_res$sample_id <- rownames(miss_svm_r_train_res); miss_svm_r_train_res <- miss_svm_r_train_res[,-c(1)]
## svm_p
miss_svm_p_train_res$sample_id <- rownames(miss_svm_p_train_res); miss_svm_p_train_res <- miss_svm_p_train_res[,-c(1)]
## pen_reg
miss_pen_reg_train_res$sample_id <- rownames(miss_pen_reg_train_res); miss_pen_reg_train_res <- miss_pen_reg_train_res[,-c(1)]
## knn
miss_knn_train_res$sample_id <- rownames(miss_knn_train_res); miss_knn_train_res <- miss_knn_train_res[,-c(1)]
## xgb
miss_xgb_train_res$sample_id <- rownames(miss_xgb_train_res); miss_xgb_train_res <- miss_xgb_train_res[,-c(1)]

# Missclasified all train. Convert into df and add method info
## rf
miss_rf_all_train<- data.frame(sample_id = miss_rf_all_train);miss_rf_all_train$RF <- rep(TRUE,nrow(miss_rf_all_train))
## svm_R
miss_svm_r_all_train <- data.frame(sample_id = miss_svm_r_all_train);miss_svm_r_all_train$SVM_R <- rep(TRUE,nrow(miss_svm_r_all_train))
## svm_p
miss_svm_p_all_train <- data.frame(sample_id = miss_svm_p_all_train);miss_svm_p_all_train$SVM_P <- rep(TRUE,nrow(miss_svm_p_all_train))
## pen_reg
miss_pen_reg_all_train <- data.frame(sample_id = miss_pen_reg_all_train);miss_pen_reg_all_train$PEN_REG <- rep(TRUE,nrow(miss_pen_reg_all_train))
## knn
miss_knn_all_train <- data.frame(sample_id = miss_knn_all_train);miss_knn_all_train$KNN <- rep(TRUE,nrow(miss_knn_all_train))
## xgb
miss_xgb_all_train <- data.frame(sample_id = miss_xgb_all_train);miss_xgb_all_train$XGB <- rep(TRUE,nrow(miss_xgb_all_train))

# Missclasified test. Convert into df and add method info
## rf
miss_rf_test<- data.frame(sample_id = miss_rf_test);miss_rf_test$RF <- rep(TRUE,nrow(miss_rf_test))
## svm_r
miss_svm_r_test <- data.frame(sample_id = miss_svm_r_test);miss_svm_r_test$SVM_R <- rep(TRUE,nrow(miss_svm_r_test))
## svm_p
miss_svm_p_test <- data.frame(sample_id = miss_svm_p_test);miss_svm_p_test$SVM_P <- rep(TRUE,nrow(miss_svm_p_test))
## pen_reg
miss_pen_reg_test <- data.frame(sample_id = miss_pen_reg_test);miss_pen_reg_test$PEN_REG <- rep(TRUE,nrow(miss_pen_reg_test))
## knn
miss_knn_test <- data.frame(sample_id = miss_knn_test);miss_knn_test$KNN <- rep(TRUE,nrow(miss_knn_test))
## xgb
miss_xgb_test <- data.frame(sample_id = miss_xgb_test);miss_xgb_test$XGB <- rep(TRUE,nrow(miss_xgb_test))

# Missclasified visualization
# Train
list_train <- list(RF = miss_rf_train_res, SVM_R = miss_svm_r_train_res, SVM_p=miss_svm_p_train_res,
                   PEN_REG=miss_pen_reg_train_res, KNN=miss_knn_train_res, XGB=miss_xgb_train_res)
list_all_train <- list(RF = miss_rf_all_train, SVM_R = miss_svm_r_all_train, SVM_p=miss_svm_p_all_train,
                       PEN_REG=miss_pen_reg_all_train, KNN=miss_knn_all_train, XGB=miss_xgb_all_train)
list_test <- list(RF = miss_rf_test, SVM_R = miss_svm_r_test, SVM_p=miss_svm_p_test,
                       PEN_REG=miss_pen_reg_test, KNN=miss_knn_test, XGB=miss_xgb_test)

for (name in c("train", "all_train", "test")) {
  lst <- get(paste("list_",name,sep=""))
  miss_train <- lst %>% reduce(full_join, by="sample_id")
  copd_ctrl_condition <- subset(combat_edata_ph, sample_id %in%
                                  miss_train$sample_id)[,pheno_cctics]
  miss_train <- merge(miss_train, copd_ctrl_condition, by="sample_id")
  if(name == "train"){
    miss_train <- miss_train %>% rename(
      Freq_RF = Freq.x,
      Freq_SVM_r = Freq.y,
      Freq_SVM_p = Freq.x.x,
      Freq_pen_reg = Freq.y.y,
      Freq_knn = Freq.x.x.x,
      Freq_xgb = Freq.y.y.y,
    )
  }
  
  # Save miss_train table
  write.csv(miss_train, paste("../../data/ML_out/two_class/",methodology,"miss_samples/miss_",name,".csv", sep=""))
  
  # Save intersection samples for further analysis
  assign(paste("miss_intersection_",name,sep=""),miss_train[complete.cases(miss_train[, 2:7]), ])
  
  miss_train_plot <- miss_train[,c(2:7)]
  miss_train_plot[is.na(miss_train_plot)] <- 0
  
  #Save order or pheatmap clustering
  pheatm <- pheatmap(as.matrix(miss_train_plot), cluster_cols = FALSE)
  miss_train_plot$sample_id <- miss_train$sample_id  
  miss_train_plot <- merge(miss_train_plot, copd_ctrl_condition, by="sample_id")
  samples_order <- miss_train_plot[pheatm$tree_row$order,]$sample_id
  
  # Reshape data from wide to long format
  df_long <- gather(miss_train_plot, key = "metric", value = "value", -sample_id, -dis_condition, -GOLD_stage, -sex, -smoker,
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
    theme_void()+
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
               scale_fill_continuous(type = "gradient",na.value="white")+
               theme(legend.position = "right", 
                     legend.direction = "horizontal",
                     legend.key.size = unit(0.3, 'cm'),
                     legend.text = element_text(size=6),
                     legend.title = element_text(size=8),
                     axis.title = element_blank()) )
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
                   axis.title = element_blank()) )
    
    
  }
  plot <- heatmap + 
    side_bar_dis_condition + side_bar_GOLD_stage + side_bar_sex + side_bar_smoker + 
    side_bar_pneumocystis_colonization +  side_bar_age +
    side_bar_emphysema + side_bar_pred_fev_prebd +side_bar_pred_fev_postbd + 
    side_bar_pred_fvc_prebd + side_bar_pred_fvc_postbd + side_bar_pred_dlco +
    plot_layout(ncol = 1, heights = c(2,0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1))
  
  ggsave(paste("../../data/ML_out/two_class/",methodology,"miss_samples/histogram_",name,".pdf",sep=""),
         width = 12, height = 10)
  
  
}


# UPSET PLOT

listInput_train <- list(RF = miss_rf_train_res$sample_id, SVM_R = miss_svm_r_train_res$sample_id,
                        SVM_p=miss_svm_p_train_res$sample_id, PEN_REG=miss_pen_reg_train_res$sample_id,
                        KNN=miss_knn_train_res$sample_id, XGB=miss_xgb_train_res$sample_id)
listInput_all_train <- list(RF = miss_rf_all_train$sample_id, SVM_R = miss_svm_r_all_train$sample_id,
                            SVM_p=miss_svm_p_all_train$sample_id, PEN_REG=miss_pen_reg_all_train$sample_id,
                            KNN=miss_knn_all_train$sample_id, XGB=miss_xgb_all_train$sample_id)
listInput_test <- list(RF = miss_rf_test$sample_id, SVM_R = miss_svm_r_test$sample_id,
                            SVM_p=miss_svm_p_test$sample_id, PEN_REG=miss_pen_reg_test$sample_id,
                            KNN=miss_knn_test$sample_id, XGB=miss_xgb_test$sample_id)

for(input in c("listInput_train", "listInput_all_train", "listInput_test" )){
  input <- get(input)
  print(upset(fromList(input),nsets = 6, nintersects = NA,order.by = c("freq")))

}

# INTERSECTION SAMPLES IN PCA
load("../../data/ML_out/two_class/copd_ctrl_test.Rda")
load("../../data/ML_out/two_class/copd_ctrl_train.Rda")
samples_test <- rownames(copd_ctrl_test)
samples_train <- rownames(copd_ctrl_train)

miss_knn_test$sample_id
miss_svm_p_train_res$sample_id
knn_miss <- combat_edata_ph[miss_knn_test$sample_id,pheno_cctics]
table(knn_miss$GOLD_stage)
svm_p_miss <- combat_edata_ph[miss_svm_p_train_res$sample_id,pheno_cctics]
table(svm_p_miss$GOLD_stage)
pdata_test <- combat_edata_ph[samples_test,pheno_cctics]
table(pdata_test$GOLD_stage)
pdata_cv <- combat_edata_ph[samples_train, pheno_cctics]
table(pdata_cv$GOLD_stage)

## train + all_train + test
miss_intersection_train <- miss_intersection_train %>% rename(
  RF = Freq_RF ,
  SVM_R = Freq_SVM_r,
  SVM_P = Freq_SVM_p,
  PEN_REG = Freq_pen_reg,
  KNN = Freq_knn ,
  XGB= Freq_xgb ,
)

miss_intersection_methods <- unique(rbind(miss_intersection_train, 
                                          miss_intersection_all_train, miss_intersection_test))

# Enrichment analysis

miss_intersection_test

head(pdata)
pdata <- combat_edata_ph[samples_test,pheno_cctics]
pdata$miss_sample <- ifelse(pdata$sample_id %in% miss_intersection_test$sample_id, 1, 0)
pdata <- pdata[,-c(2,7)]

for(i in 1:(ncol(pdata)-1)){
  print(names(pdata[i]))
    if(class(pdata[[i]]) == "numeric"){
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
# head(miss_train)
# miss_train_res[miss_train_res == TRUE] <- 'train'
# miss_train_res <- rename(miss_train_res,
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
# pData(eset) <- merge(pData(eset), miss_train_res, by = "sample_id", all =TRUE)
# pData(eset)  <- merge(pData(eset), miss_all_train,  by = "sample_id", all =TRUE)
# pData(eset)  <- merge(pData(eset), miss_test, by = "sample_id", all = TRUE)
# 
# save(eset, file = "../../data/ExpressionSetObjects/eset_jon.Rda")
