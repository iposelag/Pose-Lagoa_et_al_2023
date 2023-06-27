#!/bin/Rscript
###############################################################################
############## Recopilate metrics results of ML models ########################
###############################################################################
## Recopilate metrics results in a .csv
## command: Rscript metrics.R folder_to_ml_results number_pos_class (Rscript metrics.R two_class/optm/dea/ 1)
##          number_of_positive_class: 1 copd/2 ctrl
## output: a .csv file with all the metrics and a .csv with few of them (of interest)

# Set seet (reproducibility)
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(dplyr); library(yardstick)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read commands
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide one arguments.")
} else if (length(args)==2) {
  a <- args[1]
  b <- as.numeric(args[2])
}

# a <- "two_class/optm/data_driven_expansion/"
# b <- 2 # (o 2)
# c <- "copd" # (o ctrl)
if(b == 1){
  c <- "copd"
  d <- 16
}else{
  c <- "ctrl"
  d <- 14
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("metrics"%in%list.files(paste("../../data/ML_out/",a, sep="")) == FALSE){
  dir.create(paste("../../data/ML_out/",a,"metrics", sep=""), showWarnings = FALSE)}

load(paste("../../data/ML_out/",a,"results.Rda", sep=""))

## Save different metric results
classifiers <- c("rf", "svm_r", "svm_p", "pen_reg", "knn", "xgb")
# splits <- c("train_res", "all_train", "test")
splits <- c("train_res", "test")

for(i in 1:length(splits)){

  all_metrics <- data.frame()
  split <- splits[[i]]

  for(j in 1:length(classifiers)){
    
    classif <- classifiers[j]
    metrics <- as.data.frame(results[[split]][[classif]][[b]])
    if(split != "train_res"){
      metrics <- metrics[, c(1,3,2)]
      colnames(metrics) <- c("metric", "estimate", "estimator")
      # minus_acc <- data.frame(
      #   metric = c("1-accuracy", "1-sens", "1-spec") ,
      #   estimate = c(1-metrics[1,2], 1-metrics[12,2], 1-metrics[13,2]),
      #   estimator = c(NA,NA,NA))
    }else{
      colnames(metrics) <- c("metric", "mean", "sd")
      # minus_acc <- data.frame(
      #   metric = c("1-accuracy", "1-sens", "1-spec") ,
      #   mean = c(1-metrics[1,2], 1-metrics[12,2], 1-metrics[13,2]),
      #   sd = c(NA,NA,NA))
      }

    # metrics <- rbind(metrics, minus_acc)
    metrics$classifier <- rep(classif,d)

    if(split != "train_res"){
      metrics <- metrics[, c("classifier", "metric", "estimate", "estimator")]
    }else{
      metrics <- metrics[, c("classifier", "metric", "mean", "sd")]
    }

    all_metrics <- rbind(all_metrics, metrics)
  }

  all_metrics$"%" <- all_metrics[,3]*100

  assign(paste("all_metrics_", splits[i], "_", c, sep=""), all_metrics)

  write.csv(all_metrics,
            paste("../../data/ML_out/",a,"metrics/all_metrics_",splits[i], "_", c,".csv",sep=""), 
            row.names = FALSE)
# 
#   sum_metrics <- all_metrics %>% filter( metric == "accuracy" |
#                                                 metric == "1-accuracy" |
#                                                 metric == "precision" |
#                                                 metric == "recall" |
#                                                 metric == "sens" |
#                                                 metric == "1-sens" |
#                                                 metric == "spec" |
#                                                 metric == "1-spec" |
#                                                 metric == "roc_auc" |
#                                                 metric == "mn_log_loss" )
#   assign(paste("sum_metrics_", splits[i],sep=""), all_metrics)
#   write.csv(sum_metrics,
#             paste("../../data/ML_out/",a,"metrics/sum_metrics_",splits[i],".csv",sep=""))
  
}
