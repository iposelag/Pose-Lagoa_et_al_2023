#######################################################
## SINGLE LAYER, FEED-FORWARD NEURAL NETWORK (KERAS) ##
#######################################################

## ----setup, echo=FALSE, cache=FALSE------------------------------------------------------------------------------------------------------
# Set working directory
setwd("../../data/raw/GSE47460/")
# Set seet (reproducibility)
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(keras);library(tensorflow); library(rsample); library(tidyverse);
library(themis); library(tidymodels); library(bundle)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  a = commandArgs(trailingOnly=TRUE)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

## Load Train and Test set
## Train
load("../../ML_out/two_class/copd_ctrl_train.Rda")
## Test
load("../../ML_out/two_class/copd_ctrl_test.Rda")


## Load the feature selection (list of genes)

if(a == "dea"){
  genes <- scan("../../DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  # test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
  #                               c("dis_condition",genes)]
}else if(a == "copd_genes"){
  genes <- scan("../../raw/copd_genes/copd_genes.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]
}else if(a == "union_expanssion"){
  genes <- scan("../../raw/copd_genes/copd_genes.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]

}else if(a == "omnipath_intersection"){
  genes <- scan("../../OmniPath/eset_expansion_intersection.txt", what = "character", sep = ",")
  copd_ctrl_train <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",genes)]
}else if(a == "omnipath_union"){
  genes <- scan("../../OmniPath/eset_expansion_union.txt", what = "character", sep = ",")
  copd_ctrl_train <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",genes)]
}else if(a == "data_driven"){
  genes <- scan("../../OmniPath/data_driven.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]
}

auc_loss <- metric_set(roc_auc, mn_log_loss)

# Folds
copd_ctrl_folds <- vfold_cv(train_data)

# RECIPES
## Recipe for models that required normalized data: en, knn
copd_ctrl_normalized_rec <- recipe(dis_condition ~. , data = train_data) %>%
  step_BoxCox(all_numeric_predictors()) %>%
  step_corr(all_numeric_predictors(), threshold = .85) %>%  # for correlated variables
  step_downsample(dis_condition) %>%
  step_normalize(all_predictors())
juiced_normalized <- juice(prep(copd_ctrl_normalized_rec)) # apply the pre-processing steps to datasets

# Control
ctrl_bayes <- control_bayes(no_improve = 15, verbose = TRUE, save_pred = TRUE, save_workflow = TRUE)

# MODEL SPECIFICATIONS & WORKFLOW

ini_mlp <- Sys.time()

mlp_spec <- mlp(
  hidden_units = tune(),
  epochs = tune()) %>%
  set_engine("keras") %>%
  set_mode("classification")

mlp_wf <- workflow() %>%
  add_recipe(copd_ctrl_normalized_rec) %>%
  add_model(mlp_spec)

# TUNING
parameters_mlp <- extract_parameter_set_dials(mlp_spec)
auc_acc_sens_spec <- metric_set(accuracy, yardstick::sens,  yardstick::spec,
                                precision, roc_auc, mn_log_loss)

mlp_res <- tune_bayes(
  mlp_wf,  # add workflow
  resamples = copd_ctrl_folds,
  param_info = parameters_mlp,
  initial = 2, # This could also be a grid search
  # initial = svm_r_res,
  iter = 20,
  metrics = auc_acc_sens_spec,
  control = ctrl_bayes
)

best_accuracy_mlp <- select_best(mlp_res, "accuracy") # Best based on accuracy

fin_mlp <- Sys.time()

# FINALIZE MODELS AND GET METRICS ON RESAMPLINGS AND TEST

final_mlp <- finalize_model(
  mlp_spec,
  best_accuracy_mlp
)

mlp_wf_updated <- mlp_wf %>%
  update_model(final_mlp)

mlp_fit <- mlp_wf_updated %>%
  fit(data = train_data)

# Save the final mlp model
mlp_bundle <- bundle(mlp_fit)
save(mlp_bundle, file = paste("../../ML_out/two_class/optm/",a,"/final_models/mlp_bundle.Rda", sep=""))

mlp_rs <- mlp_fit %>%
  fit_resamples(
    resamples = copd_ctrl_folds,
    metrics = auc_acc_sens_spec,
    control = ctrl_bayes # to obatin the prediction column
  )

# Save tuning times
times <- list("mlp" = fin_mlp-ini_mlp)

# RESULTS
results <- list("train_res" = NULL , "all_train" = NULL, "test" = NULL,
                "miss_samples_train_res" = NULL, "miss_samples_all_train" = NULL,
                "miss_samples_test" = NULL, "times" = times)
r_train <- list("mlp" = NULL)
r_all_train <- list("mlp" = NULL)
r_test <- list("mlp" = NULL)
s_train <- list("mlp" = NULL)
s_all_train <- list("mlp" = NULL)
s_test <- list("mlp" = NULL)

classifier <- mlp_rs
fit <- mlp_fit

# Training resamples
pred_train <- bind_rows(classifier$.predictions, .id = "column_label")

pred_train <-  pred_train %>% rename(
  fold = column_label,
  COPD = .pred_COPD,
  CTRL = .pred_CTRL,
  prediction = .pred_class,
  observation = dis_condition)
pred_train <- pred_train[,-7]

sample_identifier <- filter(pred_train,
                            observation == "COPD" & prediction == "CTRL"
                            |observation == "CTRL" & prediction == "COPD") %>% select(c(fold,.row))

missclasified_samples_train <- as.data.frame(table(sample_identifier$.row))
missclasified_samples_train <- missclasified_samples_train %>%
  rename(
    row_id = Var1)

rownames(missclasified_samples_train) <- rownames(train_data[
  as.integer(levels(missclasified_samples_train$row_id)),])

loss_auc_train <-  pred_train %>%
  group_by(fold) %>%
  auc_loss(truth = observation, COPD) %>%
  group_by(.metric) %>%
  summarise(
    mean = mean(.estimate, na.rm = TRUE),
    sd = sd(.estimate, na.rm = TRUE)
  )

all_metrics <- pred_train %>%
  group_by(fold) %>%
  conf_mat(observation, prediction) %>%
  mutate(summary_tbl = lapply(conf_mat, summary)) %>%
  unnest(summary_tbl)

metrics_train <- all_metrics %>%
  group_by(.metric) %>%
  summarise(
    mean = mean(.estimate, na.rm = TRUE),
    sd = sd(.estimate, na.rm = TRUE)
  )

metrics_train <- rbind(metrics_train, loss_auc_train)
matrix_train <- conf_mat_resampled(classifier)
r_train <- list("metrics" = metrics_train, "matrix" = matrix_train)
s_train <- missclasified_samples_train

results$train_res <- r_train
results$miss_samples_train_res <- s_train

# ALL TRAIN
pred_all_train <- fit %>%
        predict(train_data) %>%
        bind_cols(train_data)

matrix_all_train <- pred_all_train %>%
        conf_mat(truth = dis_condition, estimate = .pred_class)

metrics_all_train <- summary(pred_all_train %>%
        conf_mat(truth = dis_condition, estimate = .pred_class))

loss_auc_all_train <- fit %>%
        predict(train_data, type ="prob") %>%
        bind_cols(train_data) %>%
        auc_loss(truth = dis_condition, .pred_COPD)

metrics_all_train <- rbind(metrics_all_train, loss_auc_all_train)

misclassified_samples_all_train <- rownames(filter(pred_all_train,
                                 dis_condition == "COPD" & .pred_class == "CTRL"
                                 |dis_condition == "CTRL" & .pred_class == "COPD"))
r_all_train <- list("metrics" = metrics_all_train, "matrix" = matrix_all_train)
s_all_train <- misclassified_samples_all_train

  # TEST
pred_test <- fit %>%
        predict(copd_ctrl_test) %>%
        bind_cols(copd_ctrl_test)

matrix_test <- pred_test %>%
        conf_mat(truth = dis_condition, estimate = .pred_class)
metrics_test <- summary(pred_test %>%
        conf_mat(truth = dis_condition, estimate = .pred_class))

loss_auc_test <- fit %>%
        predict(copd_ctrl_test, type ="prob") %>%
        bind_cols(copd_ctrl_test) %>%
        auc_loss(truth = dis_condition, .pred_COPD)

metrics_test <- rbind(metrics_test, loss_auc_test)

misclassified_samples_test <- rownames(filter(pred_test,
                                 dis_condition == "COPD" & .pred_class == "CTRL"
                                 |dis_condition == "CTRL" & .pred_class == "COPD"))

r_test <- list("metrics" = metrics_test, "matrix" = matrix_test)
s_test <- misclassified_samples_test

# SAVE RESULTS ON RESULTS LIST
results$all_train <- r_all_train
results$miss_samples_all_train <- s_all_train
results$test <- r_test
results$miss_samples_test <- s_test

save(results, file=paste('../../ML_out/two_class/optm/',a,'/results_mlp.Rda', sep = ""))
