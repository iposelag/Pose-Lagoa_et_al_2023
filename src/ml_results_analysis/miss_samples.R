#!/bin/Rscript
###############################################################################
########## Recopilate missclassified samples of ML models #####################
###############################################################################
## Recopilate missclassified samples in a .csv
## command: Rscript metrics.R folder_to_ml_input_results  (Rscript metrics.R two_class/optm/dea/ )
## output: a .csv file with all the missclasified samples and a .csv with the phenotypic ctics of interest

# Set seet (reproducibility)
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(dplyr); library(yardstick); library(UpSetR); library(Biobase); library(purrr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read commands
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide one arguments.")
} else if (length(args)==1) {
  a <- args[1]
}

# a <- "two_class/optm/disgenet_expansion/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("miss_samples"%in%list.files(paste("../../data/ML_out/",a, sep="")) == FALSE){
  dir.create(paste("../../data/ML_out/",a,"miss_samples", sep=""), showWarnings = FALSE)}

load(paste("../../data/ML_out/",a,"results.Rda", sep=""))

## 
splits <- c("miss_samples_train_res", "miss_samples_test")
classifiers <- c("rf","svm_r","svm_p","pen_reg","knn", "xgb")

for(j in 1:length(splits)){
  
  split <- splits[j]
  for(i in 1:length(classifiers)){
    classifier <- classifiers[i]
    if(split  == "miss_samples_test"){
      df_miss <- data.frame(sample_id = results[[splits[j]]][[i]])
      df_miss[[classifier]] <- rep(TRUE,nrow(df_miss))
      assign(paste(split,"_",classifier, sep =""),df_miss)
    }else{
      df_miss <- results[[splits[j]]][[i]]
      df_miss$sample_id <- rownames(df_miss); df_miss <- df_miss[,-c(1)]; 
      names(df_miss)[names(df_miss) == "Freq"] <- paste0("Freq_", classifier)
      assign(paste(split,"_",classifier, sep =""),df_miss)
    }
  }
  
}

## Load phenotypic characteristics
load("../../data/ExpressionSetObjects/eset.Rda")
dim(eset)
n_samples <- ncol(eset)
n_genes <- nrow(eset)
phenotypic_ctics <- pData(eset)

interest <- c("dis_condition","GOLD_stage", "sex", "smoker", "age", 
              "pneumocystis_colonization", "emphysema", "pred_dlco", 
              "pred_fev_prebd", "pred_fev_postbd", "pred_fvc_prebd", 
              "pred_fvc_postbd", "platform_id")
phenotypic_ctics_interest <- phenotypic_ctics[, interest]
phenotypic_ctics_interest <- phenotypic_ctics_interest %>%
  mutate(sex = factor(sex),
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
phenotypic_ctics_interest$sample_id <- rownames(phenotypic_ctics_interest)

list_train <- list(RF = miss_samples_train_res_rf, SVM_R = miss_samples_train_res_svm_r, 
                   SVM_p = miss_samples_train_res_svm_p, GLM = miss_samples_train_res_pen_reg, 
                   KNN = miss_samples_train_res_knn, XGB = miss_samples_train_res_xgb)
list_test <- list(RF = miss_samples_test_rf, SVM_R = miss_samples_test_svm_r, SVM_p=miss_samples_test_svm_p,
                  GLM =miss_samples_test_pen_reg, KNN=miss_samples_test_knn, XGB=miss_samples_test_xgb)

# Join all the miss samples for all the classifiers and with the pheno data
# for cross-validation and test 
for (name in c("train", "test")) {
  
  lst <- get(paste("list_",name,sep=""))
  miss_train <- lst %>% reduce(full_join, by="sample_id")
  write.csv(miss_train, paste("../../data/ML_out/",a,"miss_samples/miss_samples_",name,".csv",sep=""), 
            row.names = FALSE)
  
  # Save all train (test) miss classified samples
  copd_ctrl_condition <- subset(phenotypic_ctics_interest, sample_id %in%
                                  miss_train$sample_id)
  miss_train_pheno <- merge(miss_train, copd_ctrl_condition, by="sample_id")
  
  # Save all train (test) miss classified samples with pheno ctics
  write.csv(miss_train_pheno, paste("../../data/ML_out/",a,"miss_samples/miss_samples_",name,"_pheno.csv",sep=""), 
            row.names = FALSE)
}


# Upset plot
listInput_train <- list(RF = miss_samples_train_res_rf$sample_id, SVM_R = miss_samples_train_res_svm_r$sample_id, 
                        SVM_p = miss_samples_train_res_svm_p$sample_id, GLM = miss_samples_train_res_pen_reg$sample_id, 
                        KNN = miss_samples_train_res_knn$sample_id, XGB = miss_samples_train_res_xgb$sample_id)
listInput_test <-  list(RF = miss_samples_test_rf$sample_id, SVM_R = miss_samples_test_svm_r$sample_id, 
                        SVM_p=miss_samples_test_svm_p$sample_id, GLM =miss_samples_test_pen_reg$sample_id, 
                        KNN=miss_samples_test_knn$sample_id, XGB=miss_samples_test_xgb$sample_id)

for(name in c("train", "test" )){
  input <- get(paste0("listInput_",name))
  pdf(file=paste("../../data/ML_out/",a,"miss_samples/upset_",name,".pdf",sep=""),
      width = 8,height = 6)
  print(upset(fromList(input),nsets = 6, nintersects = NA,order.by = c("freq")))
  dev.off()
}

glm <- miss_samples_test_pen_reg
svm_disgenet_expansion <- miss_samples_test_svm_r
svm_disgenet_expansion <- rename(svm_disgenet_expansion, svm_r_dgn_exp = "svm_r")
svm_expansion_union <- miss_samples_test_svm_r
svm_expansion_union <- rename(svm_expansion_union, svm_r_exp_union = "svm_r")
listInput_test <- list(SVM_r_dgn_exp = svm_disgenet_expansion$sample_id,
                       SVM_r_exp_union = svm_expansion_union$sample_id,
                       GLM = glm$sample_id)
input <- get(paste0("listInput_test"))
print(upset(fromList(input),nsets = 3, nintersects = NA,order.by = c("freq")))
