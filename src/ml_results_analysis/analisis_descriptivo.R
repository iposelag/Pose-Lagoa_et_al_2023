## Analysis phenotypic vbles between train and test
eset <- get(load(paste("../../data/ExpressionSetObjects/eset.Rda",sep = "" )))
load(paste("../../data/ML_out/two_class/copd_ctrl_test.Rda", sep = ""))
load(paste("../../data/ML_out/two_class/copd_ctrl_train.Rda", sep = ""))


# Required packages
library(Biobase); library(dplyr); library(ggplot2)
samples_train <- rownames(copd_ctrl_train)
samples_test <- rownames(copd_ctrl_test)

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

phenotypic_ctics_interest_train <- subset(phenotypic_ctics_interest, sample_id %in%
                                            samples_train)

phenotypic_ctics_interest_test <- subset(phenotypic_ctics_interest, sample_id %in%
                                           samples_test)

## TRAIN vs TEST
factor_vbles <- c("dis_condition","GOLD_stage", "sex", "pneumocystis_colonization", "smoker", "platform_id")
numeric_vbles <- c("age", "pred_dlco", "pred_fev_postbd", "pred_fev_prebd", "pred_fvc_postbd", "pred_fvc_prebd")
for(i in factor_vbles){
  print(i)
  # Train
  print("Train:")
  print(prop.table(table(phenotypic_ctics_interest_train[[i]]))*100)
  # Test
  print("Test:")
  print(prop.table(table(phenotypic_ctics_interest_test[[i]]))*100)
}
for(i in numeric_vbles){
  print(i)
  # Train
  print("Train ")
  print(summary(phenotypic_ctics_interest_train[[i]]))
  # Test
  print("Test:")
  print(summary(phenotypic_ctics_interest_test[[i]]))
}

## MISS VS NO MISS
methodology <-"optm/data_driven/"
# Load miss_samples table
miss_samples_all_test <- read.csv(paste("../../data/ML_out/two_class/",methodology,"miss_samples/miss_samples_test.csv", sep=""))
miss_samples_all_train <- read.csv(paste("../../data/ML_out/two_class/",methodology,"miss_samples/miss_samples_train.csv", sep=""))


# TEST
miss_samples_test <- (miss_samples_all_test %>% filter(!is.na(miss_samples_all_test$knn)))$sample_id
pheno_miss_test <- subset(phenotypic_ctics_interest_test, sample_id %in% miss_samples_test)
phenotypic_ctics_interest_test$miss_sample <- ifelse(phenotypic_ctics_interest_test$sample_id %in%
                                                       miss_samples_test, 1, 0)
print("--------------------------------------------------------------------------------")
print("TEST")
print("Factor varibles:")
for(i in factor_vbles){
  print(i)
  # Miss
  print("Missclassified samples:")
  print(prop.table(table(pheno_miss_test[[i]]))*100)
  # Test
  print("All samples:")
  print(prop.table(table(phenotypic_ctics_interest_test[[i]]))*100)
}
print("Numeric vbles:")
for(i in numeric_vbles){
  print(i)
  print("Missclassified samples:")
  print(summary(pheno_miss_test[,i]))
  print("All samples:")
  print(summary(phenotypic_ctics_interest_test[,i]))
  print(ggplot(phenotypic_ctics_interest_test, 
               aes(x=phenotypic_ctics_interest_test[,i], fill=as.factor(miss_sample))) +
          geom_density(alpha=0.9) +
          xlab(names(phenotypic_ctics_interest_test[i])) +
          ylab("Density") +
          scale_fill_manual(values=c("0" = "#C1D9D0", "1" = "#669BF0"), name="Missclasified Samples"))
}

print("--------------------------------------------------------------------------------")
print("TRAIN")
print("Factor variables:")
# TRAIN
miss_samples_train <- (miss_samples_all_train %>% filter(!is.na(miss_samples_all_train$Freq_svm_p)))$sample_id
miss_samples_train <- (miss_samples_all_train %>% filter(miss_samples_all_train$Freq_svm_p>9,na.rm = TRUE))$sample_id
pheno_miss_train <- subset(phenotypic_ctics_interest_train, sample_id %in% miss_samples_train)
phenotypic_ctics_interest_train$miss_sample <- ifelse(phenotypic_ctics_interest_train$sample_id %in%
                                                       miss_samples_train, 1, 0)
for(i in factor_vbles){
  print(i)
  # Miss
  print("Missclassified samples:")
  print(prop.table(table(pheno_miss_train[[i]]))*100)
  # Test
  print("All samples:")
  print(prop.table(table(phenotypic_ctics_interest_train[[i]]))*100)
}
print("Numeric variables:")
for(i in numeric_vbles){
  print(i)
  print("Missclassified samples:")
  print(summary(pheno_miss_train[,i]))
  print("All samples:")
  print(summary(phenotypic_ctics_interest_train[,i]))
  print(ggplot(phenotypic_ctics_interest_train, 
               aes(x=phenotypic_ctics_interest_train[,i], fill=as.factor(miss_sample))) +
          geom_density(alpha=0.9) +
          xlab(names(phenotypic_ctics_interest_train[i])) +
          ylab("Density") +
          scale_fill_manual(values=c("0" = "#C1D9D0", "1" = "#669BF0"), name="Missclasified Samples"))
}


copd_ctrl_train <- copd_ctrl_train %>% filter(rownames(copd_ctrl_train)%in% miss_samples_train)

copd_ctrl_test <- copd_ctrl_test %>% filter(rownames(copd_ctrl_test)%in% miss_samples_test)
pheno_miss_test$split <- rep("test", nrow(pheno_miss_test))
pheno_miss_train$split <- rep("train", nrow(pheno_miss_train))
proba2 <- rbind(pheno_miss_test, pheno_miss_train)
#
# annotation_col

genes <- scan("../../data/OmniPath/data_driven.txt", what = "character", sep = ",")
copd_ctrl_train <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                     c("dis_condition",genes)]
copd_ctrl_test <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                    c("dis_condition",genes)]
proba <- rbind(copd_ctrl_train, copd_ctrl_test)

## all samples
phenotypic_ctics_interest_test$split <- rep("test", nrow(phenotypic_ctics_interest_test))
phenotypic_ctics_interest_train$split <- rep("train", nrow(phenotypic_ctics_interest_train))
proba2 <- rbind(phenotypic_ctics_interest_test, phenotypic_ctics_interest_train)
correlation <- cor(t(proba[,-1]), method = "pearson")
dim(proba)

library(pheatmap)
# heatmap(correlation,
#         col = colorRampPalette(c("blue", "white", "red"))(60),
#         main = "Correlación entre muestras",cexRow = 0.5, cexCol = 0.5)
# head(correlation)
annotation_col <- data.frame(dis_condition = factor(proba2$dis_condition),
                             miss_sample = factor(proba2$miss_sample),
                             split = factor(proba2$split), GOLD_stage = factor(proba2$GOLD_stage))
                             # GOLD_stage = proba2$GOLD_stage,
                             # sex = proba2$sex,
                             # smoker =proba2$smoker,
                             # platform = proba2$platform_id)
# annotation_row <- data.frame(emphysema = as.numeric(proba2$emphysema),
#                              pred_dlco = proba2$pred_dlco,
#                              pred_fev_prebd = proba2$pred_fev_prebd,
#                              pred_fvc_prebd = proba2$pred_fvc_prebd)
rownames(annotation_col) <- proba2$sample_id
# rownames(annotation_row) <- proba2$sample_id
indices <- match(rownames(correlation), rownames(annotation_col))
annotation_col <- annotation_col[indices, , drop = FALSE]
# indices <- match(rownames(correlation), rownames(annotation_row))
# annotation_row <- annotation_row[indices, , drop = FALSE]
# ## COLORES
color_dis_condition <- c("CTRL" = "#ffdcdb", "COPD" = "#91a8d0")
color_dis_condition<- color_dis_condition[annotation_col$dis_condition]
color_GOLD_stage <- c("0-At Risk" = "#f6e0b5", "1-Mild COPD" = "#eea990",
                      "2-Moderate COPD" = "#aa6f73", "3-Severe COPD" = "#a39193",
                      "4-Very Severe COPD" = "#66545e")
color_GOLD_stage<- color_GOLD_stage[annotation_col$GOLD_stage]
color_sex <-  c("1-Male" = "#e1f7d5", "2-Female" = "#c9c9ff")
color_sex <- color_sex[annotation_col$sex]
color_smoker <- c("1-Current" = "#83adb5", "2-Ever (>100)" = "#c7bbc9",
                  "3-Never"="#5e3c58")
color_smoker <- color_smoker[annotation_col$smoker]
color_miss_sample <-  c("0" = "green","1" ="red")
color_miss_sample <- color_miss_sample[annotation_col$miss_sample]
# color_platform <- c("GPL14550" = "#ece6ff", "GPL6480" = "#efbbff")
# color_platform <- color_platform[annotation_col$platform]
pheatmap(correlation,
         annotation_col = annotation_col,
         annotation_colors = list(
           dis_condition = color_dis_condition,
           GOLD_stage = color_GOLD_stage
         ))
