rm(list = ls()) # eliminar obxectos q poida haber na memoria

## ----------------------------------------------------------------------------------------------------------------------------------------
# GRID APPROACH
## ----------------------------------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------------------------------
# FLAGS
FLAGS <- flags(
  flag_integer("dense_units_1",32), # number of neruons
  flag_integer("dense_units_2",16),
  flag_numeric("dropout_1", 0),
  flag_numeric("dropout_2", 0),
  flag_integer("batch_size", 35),
  flag_integer("epochs", 20)
)

## ----------------------------------------------------------------------------------------------------------------------------------------
# MODEL
model <- # initialize model
  keras_model_sequential() %>%
  # engadimos as capas neuronais
  layer_dense(units = FLAGS$dense_units_1,
              activation = "relu",
              input_shape = c(ncol(train_data_norm))) %>%
  layer_dropout(FLAGS$dropout_1) %>%
  layer_dense(units = FLAGS$dense_units_2,
              activation = "relu") %>%
  layer_dropout(FLAGS$dropout_2) %>%
  # layer_dense(units = 250,
  #             activation = "relu") %>%
  # a última capa do modelo ten q ter unha única neurona e unha funcion de activacion sigmoidea
  # xa que se estar realizando unha tarefa de clasificación binaria
  layer_dense(units = 1, activation = "sigmoid")


model %>% compile(
  loss = "binary_crossentropy", # binary response
  optimizer = optimizer_adam(lr = 0.001),
  metrics = c("accuracy")
)

history <- model %>% fit(
  train_data_norm, train_labels,
  epochs = FLAGS$epochs,
  batch_size = FLAGS$batch_size,
  view_metrics = TRUE,
  validation_split = 0.2,
  verbose = 1
)

results <- model %>% evaluate(
  test_data_norm, test_labels,
  verbose = 0
)
