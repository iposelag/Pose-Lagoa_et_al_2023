## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(dplyr); library(ggplot2); library(dplyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Read cv_rep or cv_simple or optm from command line
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  a = commandArgs(trailingOnly=TRUE)
}

a <- "two_class"
b <- "dea"

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories

# cv_simple <- read.csv("../../data/ML_out/two_class/cv_simple/mrmr/metrics/sum_metrics_train_res.csv")
cv_rep <- read.csv(paste0("../../data/ML_out/",a,"/cv_rep/",b,"/metrics/sum_metrics_train_res.csv"))
optm_6 <- read.csv(paste0("../../data/ML_out/",a,"/optm/",b,"/metrics/sum_metrics_train_res.csv"))
# optm_cv <- read.csv("../../data/ML_out/two_class/cv_rep/mrmr/metrics/sum_metrics_train_res.csv")

# cv_simple$methodology <- rep("cv", nrow(cv_simple)) 
cv_rep$methodology <- rep("cv_rep", nrow(cv_rep)) 
optm_6$methodology <- rep("dea", nrow(optm_6)) 
# optm_cv$methodology <- rep("bayes_cv", nrow(optm_cv)) 

# results <- rbind(cv_simple,cv_rep,optm_6, optm_cv)
results <- rbind(cv_rep,optm_6)

data_to_plot <- results %>% filter(metric =="accuracy") 

b <- "mrmr_76"

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
### TRAINS VS TEST OPTM ###
a <- "two_class"
b <- "mrmr/mrmr_256"

optm_6 <- read.csv(paste0("../../data/ML_out/",a,"/optm/",b,"/metrics/sum_metrics_train_res.csv"))
optm_test <- read.csv(paste0("../../data/ML_out/",a,"/optm/",b,"/metrics/sum_metrics_test.csv"))

optm_6$methodology <- rep("Train", nrow(optm_6)) 
optm_test$methodology <- rep("Test", nrow(optm_test)) 
# Rename columns of test set in order to join both dataframes
optm_test <- rename(optm_test, mean = estimate, sd = estimator)

# Join dataset and reshape to plot
results <- rbind(optm_6, optm_test)
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
