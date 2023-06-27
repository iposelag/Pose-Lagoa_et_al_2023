#!/bin/Rscript
###############################################################################
######## Recopilate metrics results of ML models for RadarChart plots ########
###############################################################################
## Recopilate metrics results in a .csv
## command: Rscript all_metrics_results.R positive_class (here it by default copd)
## output: a .csv file with all the metrics for all the inputs

## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(dplyr); library(ggplot2);

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("tables"%in%list.files("../../data/ML_out/two_class/optm/") == FALSE){
  dir.create("../../data/ML_out/two_class/optm/tables", showWarnings = FALSE)}

positive_class <- "copd"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Metrics

# All inputs
ml_input <- c("data_driven", "data_driven_expansion", "disgenet", "disgenet_expansion", 
              "disgenet_entire_list","expansion_intersection", "expansion_union",
              "data_driven_guildify", "disgenet_curated_guildify")

## Cross-validation
for(i in 1:length(ml_input)){
        results_train <- read.csv(
          paste0("../../data/ML_out/two_class/optm/",ml_input[i],
                 "/metrics/all_metrics_train_res_",positive_class,".csv"))
        results_train$positive_class <- rep(positive_class, nrow(results_train))
        results_train$ml_input <- rep(ml_input[i], nrow(results_train))
        
        ## Add normalized mcc
        mcc_rows <- subset(results_train, metric == 'mcc')
        # Compute the normalized MCC
        mcc_rows$mean <- (mcc_rows$mean + 1) / 2
        mcc_rows$metric <- "mcc_normalized"
        # Combine the new row with the original data frame
        results_train <- rbind(results_train, mcc_rows)
        
        assign(paste(ml_input[i],"_results_train", sep=""), results_train)
}

results_train<- rbind(disgenet_results_train, disgenet_expansion_results_train,
                      disgenet_entire_list_results_train,
                      data_driven_results_train, data_driven_expansion_results_train,
                      expansion_intersection_results_train, expansion_union_results_train,
                      data_driven_guildify_results_train, disgenet_curated_guildify_results_train)

write.csv(results_train, file = "../../data/ML_out/two_class/optm/tables/results_train.csv", row.names = FALSE)

## Test
for(i in 1:length(ml_input)){
  results_test <- read.csv(
    paste0("../../data/ML_out/two_class/optm/",ml_input[i],
           "/metrics/all_metrics_test_",positive_class,".csv"))
  results_test$positive_class <- rep(positive_class, nrow(results_test))
  results_test$ml_input <- rep(ml_input[i], nrow(results_test))
  
  ## Add normalized mcc
  mcc_rows <- subset(results_test, metric == 'mcc')
  # Compute the normalized MCC
  mcc_rows$estimate <- (mcc_rows$estimate + 1) / 2
  mcc_rows$metric <- "mcc_normalized"
  # Combine the new row with the original data frame
  results_test <- rbind(results_test, mcc_rows)
  assign(paste(ml_input[i],"_results_test", sep=""), results_test)
}

results_test<- rbind(disgenet_results_test, disgenet_expansion_results_test,
                     disgenet_entire_list_results_test,
                      data_driven_results_test, data_driven_expansion_results_test,
                      expansion_intersection_results_test, expansion_union_results_test,
                     data_driven_guildify_results_test, disgenet_curated_guildify_results_test)

write.csv(results_test, file = "../../data/ML_out/two_class/optm/tables/results_test.csv", row.names = FALSE)

