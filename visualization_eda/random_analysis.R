###########################################
## analysis of random ##
###########################################
## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(sva); library(factoextra); library(dplyr); library(Biobase);

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data
all_results_train <- data.frame()
all_results_test <- data.frame()
all_results_time <- data.frame()
labels <- data.frame(classifier = c("RF", "SVM-rad", "SVM-poly", "GLM", "kNN", "XGB"))
# labels <- data.frame(classifier = c("rf"))

# set the folder path
folder_path <- "../../data/ML_out/two_class/optm/random_1553//"

# get the file names that match the pattern
file_names <- list.files(folder_path, pattern = "results_\\d+\\.Rda")

# read in the files using lapply and assign them to a list

# for( i in seq(10, 500, by=5)){
for( i in 1:length(file_names)){
  
  # labels$it <- rep(i,6)
  load(paste0(folder_path,file_names[i]))
  
  accuracy_train <- results$train_res
  accuracy_test <- results$test
  time_tunning <- as.data.frame(results$times)
  iteration <- data.frame(iteration  = i)
  
  classifiers_train <- rbind(
    accuracy_train$rf$metrics %>% filter(.metric == "accuracy"),
    accuracy_train$svm_r$metrics %>% filter(.metric == "accuracy"),
    accuracy_train$svm_p$metrics %>% filter(.metric == "accuracy"),
    accuracy_train$pen_reg$metrics %>% filter(.metric == "accuracy"),
    accuracy_train$knn$metrics %>% filter(.metric == "accuracy"),
    accuracy_train$xgb$metrics %>% filter(.metric == "accuracy")
  )
  
  classifiers_test <- rbind(
    accuracy_test$rf$metrics %>% filter(.metric == "accuracy"),
    accuracy_test$svm_r$metrics %>% filter(.metric == "accuracy"),
    accuracy_test$svm_p$metrics %>% filter(.metric == "accuracy"),
    accuracy_test$pen_reg$metrics %>% filter(.metric == "accuracy"),
    accuracy_test$knn$metrics %>% filter(.metric == "accuracy"),
    accuracy_test$xgb$metrics %>% filter(.metric == "accuracy")
  )
  
  classifiers_time <- cbind(iteration, time_tunning)
  all_results_time <- rbind(all_results_time, classifiers_time)
  
  classifiers_train <- cbind(labels, classifiers_train)
  all_results_train <- rbind(all_results_train, classifiers_train)
  
  classifiers_test <- cbind(labels, classifiers_test)
  all_results_test <- rbind(all_results_test, classifiers_test)
  
}

## DisGeNet
# results_disgenet <- get(load("../../data/ML_out/two_class/optm/disgenet/results.Rda"))
# rf_acc <- results_disgenet$train_res$rf$metrics %>% filter(.metric == "accuracy")
# knn_acc <- results_disgenet$train_res$knn$metrics %>% filter(.metric == "accuracy")
# pen_reg_acc <- results_disgenet$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
# svm_p_acc <- results_disgenet$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
# svm_r_acc <- results_disgenet$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
# xgb_acc <- results_disgenet$train_res$xgb$metrics %>% filter(.metric == "accuracy")

## DEA
# results_dea <- get(load("../../data/ML_out/two_class/optm/dea/results.Rda"))
# rf_acc <- results_dea$train_res$rf$metrics %>% filter(.metric == "accuracy")
# knn_acc <- results_dea$train_res$knn$metrics %>% filter(.metric == "accuracy")
# pen_reg_acc <- results_dea$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
# svm_p_acc <- results_dea$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
# svm_r_acc <- results_dea$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
# xgb_acc <- results_dea$train_res$xgb$metrics %>% filter(.metric == "accuracy")

# ## mrmr_100
# results_mrmr <- get(load("../../data/ML_out/two_class/optm/mrmr/mrmr_100/results.Rda"))
# rf_acc <- results_mrmr$train_res$rf$metrics %>% filter(.metric == "accuracy")
# knn_acc <- results_mrmr$train_res$knn$metrics %>% filter(.metric == "accuracy")
# pen_reg_acc <- results_mrmr$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
# svm_p_acc <- results_mrmr$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
# svm_r_acc <- results_mrmr$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
# xgb_acc <- results_mrmr$train_res$xgb$metrics %>% filter(.metric == "accuracy")

# ## Data_Driven
# results_data_driven <- get(load("../../data/ML_out/two_class/optm/data_driven/results.Rda"))
# rf_acc <- results_data_driven$train_res$rf$metrics %>% filter(.metric == "accuracy")
# knn_acc <- results_data_driven$train_res$knn$metrics %>% filter(.metric == "accuracy")
# pen_reg_acc <- results_data_driven$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
# svm_p_acc <- results_data_driven$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
# svm_r_acc <- results_data_driven$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
# xgb_acc <- results_data_driven$train_res$xgb$metrics %>% filter(.metric == "accuracy")

# ## Omnipath_intersection
# results_omnipath_intersection <- get(load("../../data/ML_out/two_class/optm/omnipath_intersection/results.Rda"))
# rf_acc <- results_omnipath_intersection$train_res$rf$metrics %>% filter(.metric == "accuracy")
# knn_acc <- results_omnipath_intersection$train_res$knn$metrics %>% filter(.metric == "accuracy")
# pen_reg_acc <- results_omnipath_intersection$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
# svm_p_acc <- results_omnipath_intersection$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
# svm_r_acc <- results_omnipath_intersection$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
# xgb_acc <- results_omnipath_intersection$train_res$xgb$metrics %>% filter(.metric == "accuracy")

# ## Dgn expansion
results_dgn_expansion <- get(load("../../data/ML_out/two_class/optm/dgn_expansion/results.Rda"))
rf_acc <- results_dgn_expansion$train_res$rf$metrics %>% filter(.metric == "accuracy")
knn_acc <- results_dgn_expansion$train_res$knn$metrics %>% filter(.metric == "accuracy")
pen_reg_acc <- results_dgn_expansion$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
svm_p_acc <- results_dgn_expansion$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
svm_r_acc <- results_dgn_expansion$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
xgb_acc <- results_dgn_expansion$train_res$xgb$metrics %>% filter(.metric == "accuracy")

## UNION expansion
results_union_expansion <- get(load("../../data/ML_out/two_class/optm/omnipath_union/results.Rda"))
rf_acc <- results_union_expansion$train_res$rf$metrics %>% filter(.metric == "accuracy")
knn_acc <- results_union_expansion$train_res$knn$metrics %>% filter(.metric == "accuracy")
pen_reg_acc <- results_union_expansion$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
svm_p_acc <- results_union_expansion$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
svm_r_acc <- results_union_expansion$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
xgb_acc <- results_union_expansion$train_res$xgb$metrics %>% filter(.metric == "accuracy")


## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot distribution values

# Create a histogram
ggplot(all_results_train, aes(x = mean)) +
  geom_histogram(binwidth = 0.01, fill = "gray", color = "black") +
  facet_wrap(vars(classifier))+
  theme_minimal()+
  geom_vline(data = subset(all_results_train, classifier == "RF"), 
             aes(xintercept = rf_acc$mean), color = "#6b3e26", linewidth = 0.6) +
  geom_vline(data = subset(all_results_train, classifier == "kNN"), 
             aes(xintercept = knn_acc$mean), color = "#ffa31a", linewidth = 0.6) +  
  geom_vline(data = subset(all_results_train, classifier == "GLM"), 
             aes(xintercept = pen_reg_acc$mean), color = "#ffdc73", linewidth = 0.7) +
  geom_vline(data = subset(all_results_train, classifier == "SVM-poly"), 
             aes(xintercept = svm_p_acc$mean), color = "#bdd1a0", linewidth = 0.6) +
  geom_vline(data = subset(all_results_train, classifier == "SVM-rad"), 
             aes(xintercept = svm_r_acc$mean), color = "#ffaed7", linewidth = 0.6) +
  geom_vline(data = subset(all_results_train, classifier == "XGB"), 
             aes(xintercept = xgb_acc$mean), color = "#ff6f69", linewidth = 0.6) +
  ggtitle("Distribution of Accuracy. Union expansion results in colors.") +
  xlab("Accuracy distributions") +
  ylab("")

ggsave("../../data/ML_out/two_class/optm/plots/random/omnipath_union.pdf", width = 15, height = 8)


## Obtain a p value:()
acc_real <- c("rf_acc", "svm_r_acc", "svm_p_acc", "pen_reg_acc", "knn_acc", "xgb_acc")
acc_real <- c(rf_acc$mean, svm_r_acc$mean, svm_p_acc$mean, pen_reg_acc$mean, knn_acc$mean, xgb_acc$mean)
j <- 0
for(i in c("RF", "SVM-rad", "SVM-poly", "GLM", "kNN", "XGB")){
  j <- j+1
  print(i)
  random <- all_results_train %>% filter(classifier == i)
  random_mean <- mean(random$mean)
  print(random_mean)
  print(t.test( random$mean, mu = acc_real[j]), alternative = "greater")
  
}

random_mean <- mean(random_accuracies)
