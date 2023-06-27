#!/bin/Rscript
###############################################################################
################## Tuning methodologies comparison ############################
###############################################################################
## This script has two parts:
## 1. A comparison between tuning methodologies (grid search and cross validation)
## 2- A comparison between train and test of Bayes tuning methodology
## command: Rscript tuning_methodologies_comparison.R folder_to_ml_results
## (Rscript tuning_methodologies_comparison.R two_class dea)
## output: a .csv file with all the metrics and a .csv with few of them (of interest)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(dplyr); library(ggplot2); library(dplyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read gs_rep or gs_simple or optm from command line
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  a = commandArgs(trailingOnly=TRUE)
}

# a <- "two_class"
# b <- "dea"

## ----------------------------------------------------------------------------------------------------------------------------------------
# BAYES vs GRID SEARCH

# gs_simple <- read.csv("../../data/ML_out/two_class/gs_simple/mrmr/metrics/sum_metrics_train_res.csv")
gs_rep <- read.csv(paste0("../../data/ML_out/",a,"/gs_rep/",b,"/metrics/sum_metrics_train_res.csv"))
optm_6 <- read.csv(paste0("../../data/ML_out/",a,"/optm/",b,"/metrics/sum_metrics_train_res.csv"))
# optm_gs <- read.csv("../../data/ML_out/two_class/gs_rep/mrmr/metrics/sum_metrics_train_res.csv")

# gs_simple$methodology <- rep("cv", nrow(gs_simple))
gs_rep$methodology <- rep("gs_rep", nrow(gs_rep))
optm_6$methodology <- rep("bayes", nrow(optm_6))
# optm_gs$methodology <- rep("bayes_gs", nrow(optm_gs))

# results <- rbind(gs_simple,gs_rep,optm_6, optm_gs)
results <- rbind(gs_rep,optm_6)

data_to_plot <- results %>% filter(metric =="accuracy")

# b <- "mrmr_76"

ggplot(data_to_plot) +
  aes(
    x = methodology,
    fill = classifier,
    group = classifier,
    weight = mean
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
    title = paste(a,b, sep = " "),
    subtitle = "Train"
  ) +
  geom_text(
    aes(label = round(mean, 3), y =mean),
    hjust = -0.5,
    size = 2,
position = position_dodge(width = 0.9)
  ) +
  coord_flip() +
  theme_minimal() +
  ylim(0,1)

ggsave(paste("../../data/ML_out/",a,"/plots/metrics_comparison/",b,".pdf",sep=""),
       width = 12, height = 10)

## ----------------------------------------------------------------------------------------------------------------------------------------
### TRAINS VS TEST bayes ###
# a <- "two_class"
# b <- "mrmr/mrmr_256"

bayes_train <- read.csv(paste0("../../data/ML_out/",a,"/optm/",b,"/metrics/sum_metrics_train_res.csv"))
bayes_test <- read.csv(paste0("../../data/ML_out/",a,"/optm/",b,"/metrics/sum_metrics_test.csv"))

bayes_train$methodology <- rep("Train", nrow(bayes_train))
bayes_test$methodology <- rep("Test", nrow(bayes_test))

# Rename columns of test set in order to join both dataframes
bayes_test <- rename(bayes_test, mean = estimate, sd = estimator)

# Join dataset and reshape to plot
results <- rbind(bayes_train, bayes_test)
data_to_plot <- results %>% filter(metric =="accuracy")

b <- "mrmr_100"
ggplot(data_to_plot) +
  aes(
    x = methodology,
    fill = classifier,
    group = classifier,
    weight = mean
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
    title = "Train vs Test bayes optimization",
    subtitle = b
  ) +
  geom_text(
    aes(label = round(mean, 3), y =mean),
    hjust = -0.5,
    size = 2,
    position = position_dodge(width = 0.9)
  ) +
  coord_flip() +
  theme_minimal() +
  ylim(0,1)

ggsave(paste("../../data/ML_out/",a,"/optm/plots/metrics_comparison/",b,".pdf",sep=""),
       width = 12, height = 10)
