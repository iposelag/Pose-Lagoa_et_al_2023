#!/bin/Rscript
###############################################################################
############################ Random analysis  ################################
###############################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read gs_rep or gs_simple or optm from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide one arguments.")
} else if (length(args)==1) {
  random_num_genes <- args[1]
}

random_num_genes <- "random_1553/"

if(random_num_genes == "random_30/"){
  input_genes <- "disgenet"
}else if(random_num_genes == "random_76/"){
  input_genes <- "dea"
}else if(random_num_genes == "random_100/"){
  input_genes <- "mrmr/mrmr_100"
}else if(random_num_genes == "random_163/"){
  input_genes <- "data_driven"
}else if(random_num_genes == "random_256/"){
  input_genes <- "expansion_intersection"
}else if(random_num_genes == "random_744/"){
  input_genes <- "data_driven_expansion"
}else if(random_num_genes == "random_1065/"){
  input_genes <- "disgenet_expansion"
}else (random_num_genes == "random_1553/"){
  input_genes <- "expansion_union"
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(sva); library(factoextra); library(dplyr); library(Biobase);

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

## Random results
all_results_train <- data.frame()
all_results_test <- data.frame()
all_results_time <- data.frame()
labels <- data.frame(classifier = c("RF", "SVM-rad", "SVM-poly", "GLM", "kNN", "XGB"))

# set the folder path
folder_path <- paste0("../../data/ML_out/two_class/optm/",random_num_genes)

# get the file names that match the pattern
file_names <- list.files(folder_path, pattern = "results_\\d+\\.Rda")

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

## Real results
results <- get(load(paste0("../../data/ML_out/two_class/optm/",input_genes,"/results.Rda")))
rf_acc <- results$train_res$rf$metrics %>% filter(.metric == "accuracy")
knn_acc <- results$train_res$knn$metrics %>% filter(.metric == "accuracy")
pen_reg_acc <- results$train_res$pen_reg$metrics %>% filter(.metric == "accuracy")
svm_p_acc <- results$train_res$svm_p$metrics %>% filter(.metric == "accuracy")
svm_r_acc <- results$train_res$svm_r$metrics %>% filter(.metric == "accuracy")
xgb_acc <- results$train_res$xgb$metrics %>% filter(.metric == "accuracy")


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
  ggtitle("Distribution of Accuracy. Real results in colors.") +
  xlab("Accuracy distributions") +
  ylab("")

ggsave(paste0("../../data/ML_out/two_class/optm/plots/random/",input_genes,".pdf"), width = 15, height = 8)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Obtain a p value
