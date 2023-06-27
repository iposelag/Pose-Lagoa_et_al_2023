### ENRICHMENT
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
  mutate(dis_condition = factor(dis_condition),
         sex = factor(sex),
         pneumocystis_colonization = factor(pneumocystis_colonization),
         smoker = factor(smoker),
         GOLD_stage = factor(GOLD_stage),
         platform_id = factor(platform_id),
         age = as.numeric(age),
         emphysema = as.numeric(emphysema),
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

## MISS VS NO MISS
methodology <-"optm/data_driven/"
# Load miss_samples table
miss_samples_all_test <- read.csv(paste("../../data/ML_out/two_class/",methodology,"miss_samples/miss_samples_test.csv", sep=""))
miss_samples_all_train <- read.csv(paste("../../data/ML_out/two_class/",methodology,"miss_samples/miss_samples_train.csv", sep=""))

# TEST
miss_samples_test <- (miss_samples_all_test %>% filter(!is.na(miss_samples_all_test$svm_r)))$sample_id
pheno_miss_test <- subset(phenotypic_ctics_interest_test, sample_id %in% miss_samples_test)
phenotypic_ctics_interest_test$miss_sample <- ifelse(phenotypic_ctics_interest_test$sample_id %in%
                                                       miss_samples_test, 1, 0)

# TRAIN
miss_samples_train <- (miss_samples_all_train %>% filter(!is.na(miss_samples_all_train$Freq_knn)))$sample_id
miss_samples_train <- (miss_samples_all_train %>% filter(miss_samples_all_train$Freq_svm_p>0,na.rm = TRUE))$sample_id
pheno_miss_train <- subset(phenotypic_ctics_interest_train, sample_id %in% miss_samples_train)
phenotypic_ctics_interest_train$miss_sample <- ifelse(phenotypic_ctics_interest_train$sample_id %in%
                                                        miss_samples_train, 1, 0)




pdata <- phenotypic_ctics_interest_train
for(i in 1:(ncol(pdata)-2)){
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

