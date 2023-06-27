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
file_names <- list.files(folder_path, pattern = "mrmr_\\d+\\/results.Rda")
file_names <- list.dirs(folder_path)
filtered_names <- file_names[!grepl("/final_models", file_names)]
filtered_names <- filtered_names[-1]

# for( i in seq(10, 500, by=5)){
for( i in 1:length(filtered_names)){

  print(filtered_names[i])
  labels$it <- rep(i,6)
  file_path <- paste0(filtered_names[i],"/results.Rda")
  if (file.exists(file_path)) {
  load(file_path)

  bal_accuracy_train <- results$train_res
  bal_accuracy_test <- results$test
  
  time_tunning <- as.data.frame(results$times)
  iteration <- data.frame(iteration  = i)

  classifiers_train <- rbind(
    bal_accuracy_train$rf$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_train$svm_r$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_train$svm_p$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_train$pen_reg$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_train$knn$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_train$xgb$metrics_copd %>% filter(.metric == "bal_accuracy")
  )

  classifiers_test <- rbind(
    bal_accuracy_test$rf$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_test$svm_r$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_test$svm_p$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_test$pen_reg$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_test$knn$metrics_copd %>% filter(.metric == "bal_accuracy"),
    bal_accuracy_test$xgb$metrics_copd %>% filter(.metric == "bal_accuracy")
  )

  classifiers_time <- cbind(iteration, time_tunning)
  all_results_time <- rbind(all_results_time, classifiers_time)

  classifiers_train <- cbind(labels, classifiers_train)
  all_results_train <- rbind(all_results_train, classifiers_train)

  classifiers_test <- cbind(labels, classifiers_test)
  all_results_test <- rbind(all_results_test, classifiers_test)
  }else{
    print(paste0("Falta", file_path))
  }
}


all_results_train %>%
  group_by(classifier) %>%
  ggplot( aes(x=it, y=mean, color = classifier)) +
  labs(x = "iteration",
       y = "accuracy",
       subtitle = "Performance on train set (repeatedcross-validation)") +
  geom_line() +
  geom_point() +
  facet_wrap( ~ classifier)


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
