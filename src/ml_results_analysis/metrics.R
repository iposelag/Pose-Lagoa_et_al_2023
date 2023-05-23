#!/bin/Rscript
###############################################################################
############## Recopilate metrics results of ML models ########################
###############################################################################
## Recopilate metrics results in a .csv
## command: Rscript metrics.R folder_to_ml_results (Rscript metrics.R two_class/optm/dea/)
## output: a .csv file with all the metrics and a .csv with few of them (of interest)

# Set seet (reproducibility)
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(dplyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read cv_rep or cv_simple or optm from command line
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  a = commandArgs(trailingOnly=TRUE)
}

# a <- "two_class/optm/dea/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("metrics"%in%list.files(paste("../../data/ML_out/",a, sep="")) == FALSE){
  dir.create(paste("../../data/ML_out/",a,"metrics", sep=""), showWarnings = FALSE)}

load(paste("../../data/ML_out/",a,"results.Rda", sep=""))

## Save different metric results
classifiers <- c("rf", "svm_r", "svm_p", "pen_reg", "knn", "xgb")
splits <- c("train_res", "all_train", "test")


for(i in 1:length(splits)){

  all_metrics <- data.frame()
  split <- splits[[i]]

  for(j in 1:length(classifiers)){

    classif <- classifiers[j]
    metrics <- as.data.frame(results[[split]][[classif]][1])
    if(split != "train_res"){
      metrics <- metrics[, c(1,3,2)]
      colnames(metrics) <- c("metric", "estimate", "estimator")
      minus_acc <- data.frame(
        metric = c("1-accuracy", "1-sens", "1-spec") ,
        estimate = c(1-metrics[1,2], 1-metrics[12,2], 1-metrics[13,2]),
        estimator = c(NA,NA,NA))
    }else{
      colnames(metrics) <- c("metric", "mean", "sd")
      minus_acc <- data.frame(
        metric = c("1-accuracy", "1-sens", "1-spec") ,
        mean = c(1-metrics[1,2], 1-metrics[12,2], 1-metrics[13,2]),
        sd = c(NA,NA,NA))}

    metrics <- rbind(metrics, minus_acc)
    metrics$classifier <- rep(classif,18)

    if(split != "train_res"){
      metrics <- metrics[, c("classifier", "metric", "estimate", "estimator")]
    }else{
      metrics <- metrics[, c("classifier", "metric", "mean", "sd")]
    }

    all_metrics <- rbind(all_metrics, metrics)
  }

  all_metrics$"%" <- all_metrics[,3]*100

  assign(paste("all_metrics_", splits[i],sep=""), all_metrics)

  write.csv(all_metrics,
            paste("../../data/ML_out/",a,"metrics/all_metrics_",splits[i],".csv",sep=""))

  sum_metrics <- all_metrics %>% filter( metric == "accuracy" |
                                                metric == "1-accuracy" |
                                                metric == "precision" |
                                                metric == "recall" |
                                                metric == "sens" |
                                                metric == "1-sens" |
                                                metric == "spec" |
                                                metric == "1-spec" |
                                                metric == "roc_auc" |
                                                metric == "mn_log_loss" )
  assign(paste("sum_metrics_", splits[i],sep=""), all_metrics)
  write.csv(sum_metrics,
            paste("../../data/ML_out/",a,"metrics/sum_metrics_",splits[i],".csv",sep=""))


}
