## ----setup, echo=FALSE, cache=FALSE------------------------------------------------------------------------------------------------------
# Set working directory
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library("keras")
library(tensorflow)
library(tfruns)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read arguments from command line
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
load("../../data/ML_out/two_class/copd_ctrl_train.Rda")
## Test
load("../../data/ML_out/two_class/copd_ctrl_test.Rda")


## Load the feature selection (list of genes)
if(a == "dea"){
  genes <- scan("../../data/DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]
}else if(a == "copd_genes"){
  genes <- scan("../../data/raw/copd_genes/copd_genes.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]
}else if(a == "data_driven_expansion"){
  genes <- scan("../../data/OmniPath/expansion_data_driven_genes.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]

}else if(a == "omnipath_intersection"){
  genes <- scan("../../data/OmniPath/eset_expansion_intersection.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                       c("dis_condition",genes)]
}else if(a == "omnipath_union"){
  genes <- scan("../../data/OmniPath/eset_expansion_union.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                       c("dis_condition",genes)]
}else if(a == "data_driven"){
  genes <- scan("../../data/OmniPath/data_driven.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                  c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                                c("dis_condition",genes)]
}else if(a == "disgenet"){
  genes <- scan("../../data/raw/copd_genes/COPDDisgenet.txt", what = "character", sep = ",")
  train_data <- copd_ctrl_train[,colnames(copd_ctrl_train) %in%
                                       c("dis_condition",genes)]
  test_data <- copd_ctrl_test[,colnames(copd_ctrl_test) %in%
                               c("dis_condition",genes)]
}else if(a == "dgn_expansion"){
  genes <- scan("../../data/OmniPath/expansion_dgn_genes.txt", what = "character", sep = ",")
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

# ## ----------------------------------------------------------------------------------------------------------------------------------------
# # MLP model
model <- # initialize model
   keras_model_sequential() %>%
   # engadimos as capas neuronais
   layer_dense(units = 32, activation = "relu", input_shape = c(ncol(train_data_norm))) %>%
   # layer_dropout(0.1) %>%
   layer_dense(units = 16, activation = "relu") %>%
   # layer_dense(units = 1250, activation = "relu") %>%
   # layer_dense(units = 625, activation = "relu") %>%

#   # a última capa do modelo ten q ter unha única neurona e unha funcion de activacion sigmoidea
#   # xa que se estar realizando unha tarefa de clasificación binaria
  layer_dense(units = 1, activation = "sigmoid")
#
# model <- # initialize model
#   keras_model_sequential() %>%
#   # engadimos as capas neuronais
#   layer_dense(units = 1000, activation = "relu", input_shape = c(ncol(train_data_norm))) %>%
#   layer_dense(units = 500, activation = "relu") %>%
#   layer_dense(units = 250, activation = "relu") %>%
#   # layer_dense(units = 625, activation = "relu") %>%
#   # layer_dropout(0.2) %>%
#   # a última capa do modelo ten q ter unha única neurona e unha funcion de activacion sigmoidea
#   # xa que se estar realizando unha tarefa de clasificación binaria
#   layer_dense(units = 1, activation = "sigmoid")
#
model %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(lr = 0.001),
  metrics = c("accuracy")
)
#
# # En este ejemplo, he compilado el modelo usando una función de pérdida binaria,
# # un optimizador Adam con una tasa de aprendizaje de 0.001 y la métrica de precisión
# # (accuracy) para evaluar el rendimiento del modelo.
#
history <- model %>% fit(
  train_data_norm, train_labels,
  epochs = 200,
  batch_size = 70,
  validation_split = 0.1,
  verbose = 1
)
save(history, file = paste("../../data/ML_out/two_class/optm/",a,"/final_models/results_mlp.Rda", sep=""))
plot(history)
#
#
results <- model %>% evaluate(test_data_norm, test_labels)
print("Results test:")
print(results)
#
# input_data <- layer_input(shape = c(16235))


###############
# runs <- tuning_run("experiment.R",
#                     flags = list(dense_units_1 = c(32, 64),
#                                  dense_units_2 = c(16, 32),
#                                  dropout_1 = c(0, 0.1),
#                                  dropout_2  = c(0, 0.1),
#                                  batch_size = c(35, 70),
#                                  epochs = c(20, 100, 250)
#                                  )
#                     )
# head(runs)
