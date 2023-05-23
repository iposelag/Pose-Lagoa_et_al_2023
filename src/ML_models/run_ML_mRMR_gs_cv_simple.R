## ----setup, echo=FALSE, cache=FALSE------------------------------------------------------------------------------------------------------
# Set working directory
setwd("../../data/raw/GSE47460/")
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

## --------------------------------------------------------------------------------------------------------------------
# Create required directories
if ("ML_out"%in%list.files("../../") == FALSE){dir.create("../../ML_out", showWarnings = FALSE)}
if ("two_class"%in%list.files("../../ML_out/") == FALSE){
  dir.create("../../ML_out/two_class", showWarnings = FALSE)}
if ("cv_simple"%in%list.files("../../ML_out/two_class") == FALSE){
  dir.create("../../ML_out/two_class/cv_simple", showWarnings = FALSE)}
if ("copd_genes"%in%list.files("../../ML_out/two_class/cv_simple") == FALSE){
  dir.create("../../ML_out/two_class/cv_simple/copd_genes", showWarnings = FALSE)}
if ("res_objects"%in%list.files("../../ML_out/two_class/cv_simple/") == FALSE){
  dir.create("../../ML_out/two_class/cv_simple/res_objects", showWarnings = FALSE)}
if ("mrmr_dea"%in%list.files("../../ML_out/two_class/cv_simple/") == FALSE){
    dir.create("../../ML_out/two_class/cv_simple/mrmr_dea", showWarnings = FALSE)}
if ("mrmr"%in%list.files("../../ML_out/two_class/cv_simple/") == FALSE){
      dir.create("../../ML_out/two_class/cv_simple/mrmr", showWarnings = FALSE)}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read arguments from commnad line
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  a = commandArgs(trailingOnly=TRUE)
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Upload data
ini <- Sys.time()
load("../../ExpressionSetObjects/eset.Rda")

# BATCH effect correction (platform): ComBAT
batch = pData(eset)$platform_id
modcombat = model.matrix(~ dis_condition, data=pData(eset))
combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)

# Data frame of interest with ComBat data
combat_edata_ph <- as.data.frame(t(combat_edata))
dim(combat_edata_ph)

# Add disease condition (our classif task)
combat_edata_ph$dis_condition <- as.factor(pData(eset)$dis_condition)
combat_copd_ctrl <- combat_edata_ph
combat_copd_ctrl$dis_condition <- relevel(combat_copd_ctrl$dis_condition,
                                          ref="COPD")
dim(combat_copd_ctrl)
combat_copd_ctrl = combat_copd_ctrl[ ,
                                     c("dis_condition",names(combat_copd_ctrl)[
                                       names(combat_copd_ctrl) != "dis_condition"])]
dim(combat_copd_ctrl)

## ----------------------------------------------------------------------------------------------------------------------------------------
# ML models
## ----------------------------------------------------------------------------------------------------------------------------------------

# Split data into train and test set
copd_ctrl_split <- initial_split(combat_copd_ctrl, strata = dis_condition)
copd_ctrl_split
# Train
copd_ctrl_train <- training(copd_ctrl_split)
# write.csv(copd_ctrl_train, "../../ML_out/two_class/mRMR/copd_ctrl_train_0.01.csv")

if(a == "deg_copd_ctrl_lfc.tsv"){
  deg_copd_ctrl <- read.table(file = '../../DEA/platform_sex/deg_copd_ctrl_lfc.tsv',
                              sep = '\t', header = TRUE)
  copd_ctrl_train <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       deg_copd_ctrl$gene_symbol]
}else if(a == "copd_genes"){
  copd_genes <- scan("../../raw/copd_genes/copd_genes.txt", what = "character", sep = ",")
  copd_ctrl_train <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",copd_genes)]

}else{
  a = as.numeric(a)
}

print("Class imbalane: train set")
round(prop.table(table(copd_ctrl_train$dis_condition))*100,3)

# Test
copd_ctrl_test <- testing(copd_ctrl_split)

# Metrics
auc_acc_sens_spec <- metric_set(accuracy, yardstick::sens,  yardstick::spec,
                                precision, roc_auc, mn_log_loss)
acc_sens_spec <- metric_set(accuracy, yardstick::sens, yardstick::spec, precision)
auc_loss <- metric_set(roc_auc, mn_log_loss)

copd_ctrl_pred <- select(copd_ctrl_train, -dis_condition)

ML_models <- function(copd_ctrl_train_mrmr){

  # Folds
  copd_ctrl_folds <- vfold_cv(copd_ctrl_train_mrmr)

  # RECIPES
  ## Recipe without normalization step
  copd_ctrl_rec <- recipe(dis_condition ~. , data = copd_ctrl_train_mrmr) %>%
    step_BoxCox(all_numeric_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = .85) %>% # for correlated variables
    step_downsample(dis_condition)  # for class imbalance
  juiced <- juice(prep(copd_ctrl_rec))

  ## Recipe for models that required normalized data: en, knn
  copd_ctrl_normalized_rec <- recipe(dis_condition ~. , data = copd_ctrl_train_mrmr) %>%
    step_BoxCox(all_numeric_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = .85) %>%  # for correlated variables
    step_downsample(dis_condition) %>%
    step_normalize(all_predictors())
  juiced_normalized <- juice(prep(copd_ctrl_normalized_rec)) # apply the pre-processing steps to datasets

  # Control needed for stacking
  ctrl_grid <- control_grid(save_pred = TRUE, save_workflow = TRUE, verbose = TRUE,)
  ctrl_res <- control_resamples(save_pred = TRUE, save_workflow = TRUE, verbose = TRUE,)

  # MODEL SPECIFICATIONS & WORKFLOW
  ## Random Forest
  rf_spec <- rand_forest(
    mtry = tune(), # num of predictor vbles to sample in each split of the tree
    trees = 1000,  # num of trees
    min_n = tune() # minimum node size
  ) %>%
    set_mode("classification") %>%
    set_engine("ranger", oob.error = TRUE, importance = "impurity",seed=1234)

  rf_wf <- workflow() %>%
    add_recipe(copd_ctrl_rec) %>%
    add_model(rf_spec)

  ## Support Vector Machines: radial
  svm_r_spec <- svm_rbf(
    cost = tune(),
    rbf_sigma = tune()
  ) %>%
    set_engine("kernlab") %>%
    set_mode("classification")

  svm_r_wf <- workflow() %>%
    add_recipe(copd_ctrl_normalized_rec) %>%
    add_model(svm_r_spec)

  ## Support Vector Machines: polynomial
  svm_p_spec <-
    svm_poly(cost = tune(), degree = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")

  svm_p_wf <- workflow() %>%
    add_recipe(copd_ctrl_normalized_rec) %>%
    add_model(svm_p_spec)

  ## Penalized Regression Model
  en_spec <- logistic_reg(penalty = tune(),
                          mixture = tune()) %>%
    set_engine("glmnet") %>%
    set_mode("classification")

  en_wf <- workflow() %>%
    add_recipe(copd_ctrl_normalized_rec) %>%
    add_model(en_spec)

  ## k-Nearest Neighbours
  knn_spec <- nearest_neighbor(
    neighbors = tune(),
    dist_power = tune(),
    weight_func = tune()) %>%
    set_engine("kknn") %>%
    set_mode("classification")

  knn_wf <- workflow() %>%
    add_recipe(copd_ctrl_normalized_rec) %>%
    add_model(knn_spec)

  ## XGBoost
  xgb_spec <- boost_tree(
      trees = 1000,
      tree_depth = tune(), min_n = tune(),
      loss_reduction = tune(),
      sample_size = tune(), mtry = tune(),
      learn_rate = tune()
    ) %>%
      set_engine("xgboost") %>%
      set_mode("classification")

    xgb_wf <- workflow() %>%
        add_recipe(copd_ctrl_rec) %>%
        add_model(xgb_spec)

  # TUNING
  ## Random Forest
  parameters_rf <- parameters(rf_spec)
  #' It is necessary to finalize the range of values of mtry() that depend on the
  #' number of columns
  parameters_rf <- rf_spec %>% parameters() %>% finalize(copd_ctrl_pred)

  rf_grid <- grid_regular(parameters_rf, levels = 5)

  rf_res <- tune_grid(
    rf_wf,  # add workflow
    resamples = copd_ctrl_folds,
    metrics = auc_acc_sens_spec,
    control = ctrl_grid,
    grid = rf_grid # grid of values to try
  )
  # save(rf_res, file="../../ML_out/two_class/cv_simple/res_objects/rf_res.Rda")

  best_loss_rf <- select_best(rf_res, "mn_log_loss") # Best based on mn_log_loss
  best_accuracy_rf <- select_best(rf_res, "accuracy") # Best based on accuracy
  best_sens_rf <- select_best(rf_res, "sens") # Best based on sensitivity
  best_spec_rf <- select_best(rf_res, "spec") # Best based on specificity
  best_auc_rf <- select_best(rf_res, "roc_auc") # Best based on auc
  best_precision_rf <- select_best(rf_res, "precision") # Best based on precision

  ## Support Vector Machines: radial
  parameters_r_svm <- parameters(svm_r_spec)
  grid_r_svm <- grid_regular(parameters_r_svm, levels = 5)

  svm_r_res <- tune_grid(
    svm_r_wf,  # add workflow
    resamples = copd_ctrl_folds,
    metrics = auc_acc_sens_spec,
    grid = grid_r_svm, # grid of values to try
    control = ctrl_grid
  )
  # save(svm_r_res, file="../../ML_out/two_class/cv_simple/res_objects/svm_r_res.Rda")

  best_loss_svm_r <- select_best(svm_r_res, "mn_log_loss") # Best based on log_loss
  best_accuracy_svm_r <- select_best(svm_r_res, "accuracy") # Best based on accuracy
  best_sens_svm_r <- select_best(svm_r_res, "sens") # Best based on sensitivity
  best_spec_svm_r <- select_best(svm_r_res, "spec") # Best based on specificity
  best_auc_svm_r <- select_best(svm_r_res, "roc_auc") # Best based on auc
  best_precision_svm_r <- select_best(svm_r_res, "precision") # Best based on precsion

  ## Support Vector Machines: polynomial
  parameters_p_svm <- parameters(svm_p_spec)
  grid_p_svm <- grid_regular(parameters_p_svm, levels = 5)

  svm_p_res <- tune_grid(
    svm_p_wf,  # add workflow
    resamples = copd_ctrl_folds,
    metrics = auc_acc_sens_spec,
    control = ctrl_grid,
    grid = grid_p_svm # grid of values to try

  )
  # save(svm_p_res, file="../../ML_out/two_class/cv_simple/res_objects/svm_p_res.Rda")

  best_loss_svm_p <- select_best(svm_p_res, "mn_log_loss") # Best based on loss
  best_accuracy_svm_p <- select_best(svm_p_res, "accuracy") # Best based on accuracy
  best_sens_svm_p <- select_best(svm_p_res, "sens") # Best based on sensitivity
  best_spec_svm_p <- select_best(svm_p_res, "spec") # Best based on specificity
  best_auc_svm_p <- select_best(svm_p_res, "roc_auc") # Best based on auc
  best_precision_svm_p <- select_best(svm_p_res, "precision") # Best based on precision

  ## Penalized Regression Model
  parameters_en <- parameters(en_spec)
  en_grid <- grid_regular(penalty(), mixture(), levels = c(mixture = 4, penalty = 10))

  en_res <- tune_grid(en_wf,
                      resamples = copd_ctrl_folds,
                      grid = en_grid,
                      control = ctrl_grid,
                      metrics = auc_acc_sens_spec
                      )
  # save(en_res, file="../../ML_out/two_class/cv_simple/res_objects/en_res.Rda")

  best_loss_en <- select_best(en_res, "mn_log_loss")   # Best based on loss
  best_accuracy_en <- select_best(en_res, "accuracy") # Best based on accuracy
  best_sens_en <- select_best(en_res, "sens") # Best based on sensitivity
  best_spec_en <- select_best(en_res, "spec") # Best based on specificity
  best_auc_en <- select_best(en_res, "roc_auc") # Best based on auc
  best_precision_en <- select_best(en_res, "precision") # Best based on precision


  ## k-Nearest Neighbours
  parameters_knn <- parameters(knn_spec)
  knn_grid <- grid_regular(parameters_knn, levels = 5)

  knn_res <- tune_grid(knn_wf,
                       resamples = copd_ctrl_folds,
                       grid = knn_grid,
                       control = ctrl_grid,
                       metrics = auc_acc_sens_spec
                       )
  # save(knn_res, file="../../ML_out/two_class/cv_simple/res_objects/knn_res.Rda")

  best_loss_knn <- select_best(knn_res, "mn_log_loss") # Best based on loss
  best_accuracy_knn <- select_best(knn_res, "accuracy") # Best based on accuracy
  best_sens_knn <- select_best(knn_res, "sens") # Best based on sensitivity
  best_spec_knn <- select_best(knn_res, "spec") # Best based on specificity
  best_auc_knn <- select_best(knn_res, "roc_auc") # Best based on auc
  best_precision_knn <- select_best(knn_res, "precision") # Best based on precision

  ## XGBoost
  xgb_grid <- grid_latin_hypercube(
    tree_depth(),
    min_n(),
    loss_reduction(),
    sample_size = sample_prop(),
    finalize(mtry(), copd_ctrl_pred),
    learn_rate(),
    size = 30
  )

  xgb_res <- tune_grid(
    xgb_wf,
    resamples = copd_ctrl_folds,
    grid = xgb_grid,
    control = ctrl_grid,
    metrics = auc_acc_sens_spec
  )

  # save(xgb_res, file="../../ML_out/two_class/cv_simple/res_objects/xgb_res.Rda")

  best_loss_xgb <- select_best(xgb_res, "mn_log_loss") # Best based on loss
  best_accuracy_xgb <- select_best(xgb_res, "accuracy") # Best based on accuracy
  best_sens_xgb <- select_best(xgb_res, "sens") # Best based on sensitivity
  best_spec_xgb <- select_best(xgb_res, "spec") # Best based on specificity
  best_auc_xgb <- select_best(xgb_res, "roc_auc") # Best based on auc
  best_precision_xgb <- select_best(xgb_res, "precision") # Best based on precision

  # FINALIZE MODELS AND GET THE LOSS VALUE
  ## Random Forest
  final_rf <- finalize_model(
    rf_spec,
    best_accuracy_rf
  )

  rf_wf_updated <- rf_wf %>%
    update_model(final_rf)

  rf_fit <- rf_wf_updated %>%
    fit(data = copd_ctrl_train_mrmr)

  rf_rs <- rf_fit %>%
    fit_resamples(
      resamples = copd_ctrl_folds,
      metrics = auc_acc_sens_spec,
      control = ctrl_res # to obatin the prediction column
    )

  ## Suppor Vector Machines: radial
  final_svm_r <- finalize_model(
    svm_r_spec,
    best_accuracy_svm_r
  )

  svm_r_wf_updated <- svm_r_wf %>%
    update_model(final_svm_r)

  svm_r_fit <- svm_r_wf_updated %>%
    fit(data = copd_ctrl_train_mrmr)

  svm_r_rs <- svm_r_fit %>%
    fit_resamples(
      resamples = copd_ctrl_folds,
      metrics = auc_acc_sens_spec,
      control = ctrl_res # to obatin the prediction column
    )

  ## Suppor Vector Machines: poly
  final_svm_p <- finalize_model(
    svm_p_spec,
    best_accuracy_svm_p
  )

  svm_p_wf_updated <- svm_p_wf %>%
    update_model(final_svm_p)

  svm_p_fit <- svm_p_wf_updated %>%
    fit(data = copd_ctrl_train_mrmr)

  svm_p_rs <- svm_p_fit %>%
    fit_resamples(
      resamples = copd_ctrl_folds,
      metrics = auc_acc_sens_spec,
      control = ctrl_res # to obatin the prediction column
    )

  ## Penalized Regression Model
  final_en <- finalize_model(
    en_spec,
    best_accuracy_en
  )

  en_wf_updated <- en_wf %>%
    update_model(final_en)

  en_fit <- en_wf_updated %>%
    fit(data = copd_ctrl_train_mrmr)

  en_rs <- en_fit %>%
    fit_resamples(
      resamples = copd_ctrl_folds,
      metrics = auc_acc_sens_spec,
      control = ctrl_res # to obatin the prediction column
    )

  ## k-Nearest Neighbours
  final_knn <- finalize_model(
    knn_spec,
    best_accuracy_knn
  )

  knn_wf_updated <- knn_wf %>%
    update_model(final_knn)

  knn_fit <- knn_wf_updated %>%
    fit(data = copd_ctrl_train_mrmr)

  knn_rs <- knn_fit %>%
    fit_resamples(
      resamples = copd_ctrl_folds,
      metrics = auc_acc_sens_spec,
      control = ctrl_res # to obatin the prediction column
    )

  ## xgboost
  final_xgb <- finalize_model(
    xgb_spec,
    best_accuracy_xgb
  )

  xgb_wf_updated <- xgb_wf %>%
    update_model(final_xgb)

  xgb_fit <- xgb_wf_updated %>%
    fit(data = copd_ctrl_train_mrmr)

  xgb_rs <- xgb_fit %>%
    fit_resamples(
      resamples = copd_ctrl_folds,
      metrics = auc_acc_sens_spec,
      control = ctrl_res # to obatin the prediction column
    )

    ## STACK model (uses boostrap)
    copd_ctrl_data_st <-
        stacks() %>%
        add_candidates(rf_res) %>%
        add_candidates(svm_r_res) %>%
        add_candidates(svm_p_res) %>%
        add_candidates(en_res) %>%
        add_candidates(knn_res) %>%
        add_candidates(xgb_res)

    copd_ctrl_model_st <-
      copd_ctrl_data_st %>%
      blend_predictions(penalty = 10^seq(-2, -0.5, length = 20)) #tune over penalty

    st_fit <-
      copd_ctrl_model_st %>%
      fit_members()

    print(st_fit)
    #  collect_parameters(copd_ctrl_model_st, "svm_r_res") %>% filter(coef != 0)


    # RESULTS
    results <- list("train_res" = NULL , "all_train" = NULL, "test" = NULL,
                    "miss_samples_train_res" = NULL, "miss_samples_all_train" = NULL,
                    "miss_samples_test" = NULL)
    r_train <- list("rf" = NULL , "svm_r" = NULL, "svm_p" = NULL, "pen_reg" = NULL,
              "knn" = NULL, "xgb" = NULL)
    r_all_train <- list("rf" = NULL , "svm_r" = NULL, "svm_p" = NULL, "pen_reg" = NULL,
              "knn" = NULL, "xgb" = NULL, "stack" = NULL)
    r_test <- list("rf" = NULL , "svm_r" = NULL, "svm_p" = NULL, "pen_reg" = NULL,
              "knn" = NULL, "xgb" = NULL, "stack" = NULL)
    s_train <- list("rf" = NULL , "svm_r" = NULL, "svm_p" = NULL, "pen_reg" = NULL,
              "knn" = NULL, "xgb" = NULL, "stack" = NULL)
    s_all_train <- list("rf" = NULL , "svm_r" = NULL, "svm_p" = NULL, "pen_reg" = NULL,
              "knn" = NULL, "xgb" = NULL, "stack" = NULL)
    s_test <- list("rf" = NULL , "svm_r" = NULL, "svm_p" = NULL, "pen_reg" = NULL,
              "knn" = NULL, "xgb" = NULL, "stack" = NULL)
    classifiers <- list(rf_rs, svm_r_rs, svm_p_rs, en_rs, knn_rs, xgb_rs)
    fits <- list(rf_fit, svm_r_fit, svm_p_fit, en_fit, knn_fit, xgb_fit, st_fit)


    # Training (aka results from resamples with better accuracy)
    for(i in 1:length(classifiers)){

            classifier <- classifiers[[i]]
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
            rownames(missclasified_samples_train) <- rownames(copd_ctrl_train_mrmr[
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
            r_train[[i]] <- list("metrics" = metrics_train, "matrix" = matrix_train)
            s_train[[i]] <- missclasified_samples_train
        }

    results$train_res <- r_train
    results$miss_samples_train_res <- s_train

    # ALL TRAIN & TEST
    for(i in 1:length(fits)){

            fit <- fits[[i]]

            # ALL TRAIN
            pred_all_train <- fit %>%
                    predict(copd_ctrl_train_mrmr) %>%
                    bind_cols(copd_ctrl_train_mrmr)

            matrix_all_train <- pred_all_train %>%
                    conf_mat(truth = dis_condition, estimate = .pred_class)

            metrics_all_train <- summary(pred_all_train %>%
                    conf_mat(truth = dis_condition, estimate = .pred_class))

            loss_auc_all_train <- fit %>%
                    predict(copd_ctrl_train_mrmr, type ="prob") %>%
                    bind_cols(copd_ctrl_train_mrmr) %>%
                    auc_loss(truth = dis_condition, .pred_COPD)

            metrics_all_train <- rbind(metrics_all_train, loss_auc_all_train)

            misclassified_samples_all_train <- rownames(filter(pred_all_train,
                                             dis_condition == "COPD" & .pred_class == "CTRL"
                                             |dis_condition == "CTRL" & .pred_class == "COPD"))

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


            r_all_train[[i]] <- list("metrics" = metrics_all_train, "matrix" = matrix_all_train)
            s_all_train[[i]] <- misclassified_samples_all_train
            r_test[[i]] <- list("metrics" = metrics_test, "matrix" = matrix_test)
            s_test[[i]] <- misclassified_samples_test


        }

    # ALL TRAIN
    results$all_train <- r_all_train
    results$miss_samples_all_train <- s_all_train
    # TEST
    results$test <- r_test
    results$miss_samples_test <- s_test


  return(results)

}

if(class(a) == "numeric"){

  for( i in seq(a, a+15, by=5)){

    if( i == 505){
      break
    }

    mrmr_copd_ctrl_train <- read.table(file = paste("../../ML_out/two_class/mRMR/mrmr/",i,"_output.txt",sep=""),
                                       sep = '\t', header = TRUE)
    mrmr_copd_ctrl_train$Name <- gsub(" ", "", mrmr_copd_ctrl_train$Name)
    copd_ctrl_train_mrmr <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                              mrmr_copd_ctrl_train$Name]
    copd_ctrl_train_mrmr$dis_condition <- copd_ctrl_train$dis_condition
    results <- ML_models(copd_ctrl_train_mrmr)

    print(paste("Results of iter", i, sep=""))
    print(results)
    # Save as a Rdata file
    save(results, file=paste("../../ML_out/two_class/cv_simple/mrmr/results_",i,".Rda", sep= ""))

  }
}else{

    results <- ML_models(copd_ctrl_train)
    print(results)
    # Save as a Rdata file
    save(results, file = paste('../../ML_out/two_class/cv_simple/',a,'/results.Rda', sep = ""))
}

fin <- Sys.time()
print(fin-ini)
