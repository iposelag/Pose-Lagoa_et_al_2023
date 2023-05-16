###########################################
## analysis of mrmr results from 10-500 ##
###########################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(sva); library(factoextra); library(dplyr); library(Biobase);

##########################3
#genes
# all_results <- data.frame()
# labels <- data.frame(classifier = c("rf", "svm_r", "svm_p", "pen_reg", "knn"))
#
#
# i=5
# for (file_name in file_list) {
#   # construct the full file path
#   i = i+5
#   labels$it <- rep(i,6)
#   file_path <- paste0("../../data/ML_out/two_class/cv_rep/random_33/", file_name)
#
#   # read the file using the appropriate function (e.g. read.csv, read.table, read_excel, etc.)
#   load(file_path)
#
#   # do something with the data (e.g. process it, analyze it, etc.)
#   accuracy_train <- results$train_res
#   accuracy_test <- results$test
#
#   classifiers_train <- rbind(
#     accuracy_train$rf$metrics %>% filter(.metric == "accuracy"),
#     accuracy_train$svm_r$metrics %>% filter(.metric == "accuracy"),
#     accuracy_train$svm_p$metrics %>% filter(.metric == "accuracy"),
#     accuracy_train$pen_reg$metrics %>% filter(.metric == "accuracy"),
#     accuracy_train$knn$metrics %>% filter(.metric == "accuracy"),
#     accuracy_train$xgb$metrics %>% filter(.metric == "accuracy")
#   )
#
#   classifiers_test <- rbind(
#     accuracy_test$rf$metrics %>% filter(.metric == "accuracy"),
#     accuracy_test$svm_r$metrics %>% filter(.metric == "accuracy"),
#     accuracy_test$svm_p$metrics %>% filter(.metric == "accuracy"),
#     accuracy_test$pen_reg$metrics %>% filter(.metric == "accuracy"),
#     accuracy_test$knn$metrics %>% filter(.metric == "accuracy"),
#     accuracy_test$xgb$metrics %>% filter(.metric == "accuracy")
#   )
#
#   classifiers_train <- cbind(labels, classifiers_train)
#   all_results_train <- rbind(all_results_train, classifiers_train)
#
#   classifiers_test <- cbind(labels, classifiers_test)
#   all_results_test <- rbind(all_results_test, classifiers_test)
# }

all_results_train <- data.frame()
all_results_test <- data.frame()
# labels <- data.frame(classifier = c("rf", "svm_r", "svm_p", "pen_reg", "knn", "xgb"))
labels <- data.frame(classifier = c("rf"))

# for( i in seq(10, 500, by=5)){
for( i in seq(2, 2000, by=2)){

  # labels$it <- rep(i,6)
  load(paste("../../data/ML_out/two_class/cv_rep/random_33/results_",i,".Rda", sep = ""))

  accuracy_train <- results$train_res
  accuracy_test <- results$test

  classifiers_train <- rbind(
    accuracy_train$rf$metrics %>% filter(.metric == "accuracy")
    # accuracy_train$svm_r$metrics %>% filter(.metric == "accuracy"),
    # accuracy_train$svm_p$metrics %>% filter(.metric == "accuracy"),
    # accuracy_train$pen_reg$metrics %>% filter(.metric == "accuracy"),
    # accuracy_train$knn$metrics %>% filter(.metric == "accuracy"),
    # accuracy_train$xgb$metrics %>% filter(.metric == "accuracy")
  )

  classifiers_test <- rbind(
    accuracy_test$rf$metrics %>% filter(.metric == "accuracy")
    # accuracy_test$svm_r$metrics %>% filter(.metric == "accuracy"),
    # accuracy_test$svm_p$metrics %>% filter(.metric == "accuracy"),
    # accuracy_test$pen_reg$metrics %>% filter(.metric == "accuracy"),
    # accuracy_test$knn$metrics %>% filter(.metric == "accuracy"),
    # accuracy_test$xgb$metrics %>% filter(.metric == "accuracy")
  )

  classifiers_train <- cbind(labels, classifiers_train)
  all_results_train <- rbind(all_results_train, classifiers_train)

  classifiers_test <- cbind(labels, classifiers_test)
  all_results_test <- rbind(all_results_test, classifiers_test)

}

load("../../data/ML_out/two_class/disgenet/results.Rda")
rf_acc <- accuracy_test$rf$metrics %>% filter(.metric == "accuracy")

# Subset the data to only include accuracy values
acc <- all_results_train[all_results_train$metric == "mean",]
hist(all_results_train$mean)

# Create a histogram
ggplot(all_results_train, aes(x = mean)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
  ggtitle("Distribution of Accuracy") +
  xlab("Accuracy") +
  ylab("Frequency")

# Or create a density plot
ggplot(acc, aes(x = mean)) +
  geom_density(fill = "lightblue", color = "black") +
  ggtitle("Distribution of Accuracy") +
  xlab("Accuracy") +
  ylab("Density")

# TRAIN

p1 <- all_results_train %>%
  group_by(classifier) %>%
  ggplot( aes(x=it, y=mean, color = classifier)) +
  scale_color_manual(
    values = c(rf = "#6b3e26",
               svm_r = "#ffc5d9",
               svm_p = "#bdd1a0",
               knn = "#ffad60",
               pen_reg = "#feda75",
               xgb = "#ff6f69"
    )) +
  labs(x = "Number of genes",
       y = "Accuracy") +
  theme(legend.position='none') +
  labs(
    title = "Two class cv_rep",
    subtitle = "Train"
  )+
  geom_line() +
  geom_point() +
  facet_wrap( ~ classifier)

p1
ggsave("../../data/ML_out/two_class/cv_rep/train.png", width = 15, height = 8)

p2 <- all_results_train %>%
  group_by(classifier) %>%
  ggplot( aes(x=it, y=mean, color = classifier)) +
  labs(x = "iteration",
       y = "accuracy",
       subtitle = "Performance on train set (repeatedcross-validation)") +
  geom_line() +
  geom_point()
p2
ggsave("../../data/Visualization/train_join_bayes.png", width = 12, height = 8)

# TEST

p1 <- all_results_test %>%
  group_by(classifier) %>%
  ggplot( aes(x=it, y=.estimate, color = classifier)) +
  labs(x = "iteration",
       y = "accuracy",
       subtitle = "Performance on test set (repeatedcross-validation)") +
  theme(legend.position='none') +
  geom_line() +
  geom_point() +
  facet_wrap( ~ classifier)+
  ylim(0.5, 1)
p1
ggsave("accuracy_grouped_test.png", width = 15, height = 8)

p2 <- all_results_test %>%
  group_by(classifier) %>%
  ggplot( aes(x=it, y=.estimate, color = classifier)) +
  labs(x = "iteration",
       y = "accuracy",
       subtitle = "Performance on test set (repeatedcross-validation)") +
  geom_line() +
  geom_point()+
  ylim(0, 1)
p2
ggsave("accuracy_test.png", width = 12, height = 8)




# erro <-results$all
#
# if(exists("erro")){
#   results$train_res <- results$train
#   results$miss_samples_train_res <- results$miss_samples_train
#   results$train <- NULL
#   results$miss_samples_train <- NULL
# }
