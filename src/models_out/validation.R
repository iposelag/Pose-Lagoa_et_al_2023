# Set seet (reproducibility)
set.seed(1234)


## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(tidyverse)
library(tidymodels)
library(Biobase)
library(vip)
library(sva)
library(themis)
library(glmnet)
library(stacks)
library(bundle)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Upload data

## Validation data
copd_ctrl_validation <- get(load("../../data/ML_out/validation_data/copd_ctrl_jon.Rda"))

## Models
input <- "data_driven"
model <- "rf_fit"
load(paste0("../../data/ML_out/two_class/optm/",input,"/final_models/",model,".Rda"))


pred_test <- model %>%
        predict(copd_ctrl_validation) %>%
        bind_cols(copd_ctrl_validation)

matrix_test <- pred_test %>%
        conf_mat(truth = dis_condition, estimate = .pred_class)
metrics_test <- summary(pred_test %>%
        conf_mat(truth = dis_condition, estimate = .pred_class))

loss_auc_test <- fit %>%
        predict(copd_ctrl_validation, type ="prob") %>%
        bind_cols(copd_ctrl_validation) %>%
        auc_loss(truth = dis_condition, .pred_COPD)

metrics_test <- rbind(metrics_test, loss_auc_test)

misclassified_samples_test <- rownames(filter(pred_test,
                                 dis_condition == "COPD" & .pred_class == "CTRL"
                                 |dis_condition == "CTRL" & .pred_class == "COPD"))

#
# r_all_train[[i]] <- list("metrics" = metrics_all_train, "matrix" = matrix_all_train)
# s_all_train[[i]] <- misclassified_samples_all_train
# r_test[[i]] <- list("metrics" = metrics_test, "matrix" = matrix_test)
# s_test[[i]] <- misclassified_samples_test