## ----setup, echo=FALSE, cache=FALSE------------------------------------------------------------------------------------------------------
# Set working directory
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library("keras")
library(tensorflow)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

## Load Train and Test set
## Train
load("../../ML_out/two_class/copd_ctrl_train.Rda")
## Test
load("../../ML_out/two_class/copd_ctrl_test.Rda")


## Load the feature selection (list of genes)
a = "dea"

if(a == "dea"){
  genes <- scan("../../DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]
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


## ----------------------------------------------------------------------------------------------------------------------------------------
# Preprocessing

## Transform the class labels into numerical labels
train_data$dis_condition <- ifelse(train_data$dis_condition == "COPD", 1, 0)
test_data$dis_condition <- ifelse(test_data$dis_condition == "COPD", 1, 0)

## Delete the class labels from the training & testing sets
## Train
train_labels <- train_data$dis_condition
train_data <- train_data[, -c(1)]
## Test
test_labels <- test_data$dis_condition
test_data <- test_data[, -c(1)]

## Noramlize training and test data
train_data_norm <- scale(train_data)
test_data_norm <- scale(test_data)

## ----------------------------------------------------------------------------------------------------------------------------------------
# BAYESIAN APPROACH
## ----------------------------------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------------------------------
##  Funtion with the model

training_credit = function(initParams){

  # Define Model
  model <- keras_model_sequential() %>%
    # capas neuronais
    layer_dense(units = initParams$dense_units_1,
      activation = 'reulu',
      input_shape = c(ncol(train_data_norm))) %>%
    layer_dropout(rate = initParams$dropout_1) %>%
    layer_dense(units = initParams$dense_units_2,
      activation = 'relu') %>%
    layer_dropout(rate = initParams$dropout_2) %>%
    layer_dense(units = 1, activation = 'sigmoid')

  summary(model)

  model %>% compile(
    loss = 'binary_crossentropy',
    optimizer = optimizer_adam(lr = 0.001),
    metrics = c('accuracy')
  )

  # Training & Evaluation

  history = model %>% fit(
    train_data_norm, train_labels,
    epochs = 30,
    batch_size = FLAGS$batch_size,
    view_metrics = TRUE,
    validation_split = 0.2,
    verbose = 1
  )

  score = model %>% evaluate(
    test_data_norm, test_labels,
    verbose = 0
  )

  A = model %>% predict_classes(x_test)

  A = as.factor(A)
  y_test2 = y_test
  y_test2 = as.factor(y_test2)

  con = confusionMatrix(A, y_test2, positive="1")

return(con$table[2,2])
}

initParams = list ( dropout_1 = 0.1, dropout_2 = 0.1)
maximizeACC = function(dropout_1, dropout_2) {

  replaceParams = list ( dropout1 = dropout_1, dropout2 = dropout_2)
  updatedParams = modifyList(initParams, replaceParams)

  score = training_credit(updatedParams)
  results = list (Score = score,  Pred = 0)
  return(results)
}


boundsParams = list(dense_units_1 = c(32, 64),
                                 dense_units_2 = c(16, 32),
                                 dropout_1 = c(0, 0.1),
                                 dropout_2  = c(0, 0.1),
                                 batch_size = c(35, 70),
                                 epochs = c(20, 100, 250)
                                 )

Final_calibrated = BayesianOptimization(maximizeACC, bounds = boundsParams,
  init_grid_dt = as.data.table(boundsParams),
    init_points = 10, n_iter = 30, acq = "ucb",
  kappa = 2.576, eps = 0, verbose = FALSE)


tail(Final_calibrated$History)
