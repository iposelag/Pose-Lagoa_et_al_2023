#!/bin/Rscript
###############################################################################
#### Barplots for the comparison among the different ML methods and using  ####
#### different inputs gene sets ###############################################
###############################################################################
## This script has two parts:
## 1. Barplot generation for the comparison between different ML methods and inputs
## for cross-validation and test sets
## 2. Barplot generation for the comparison among the number of genes by input
## 3. Dirty RadarChart
## command: Rscript barplot_ml_input_comparison.R
## output: 3 barplots (cross-validation, test comparisons and RadarChart)

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

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

## Cross-validation
dea_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/dea/metrics/all_metrics_train_res.csv"))
# minmax_100_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_100/metrics/all_metrics_train_res.csv"))
data_driven_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/data_driven/metrics/all_metrics_train_res.csv"))
# minmax_30_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_30/metrics/all_metrics_train_res.csv"))
# minmax_76_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_76/metrics/all_metrics_train_res.csv"))
# minmax_163_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_163/metrics/all_metrics_train_res.csv"))
dgn_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/disgenet/metrics/all_metrics_train_res.csv"))
expansion_intersection_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/expansion_intersection/metrics/all_metrics_train_res.csv"))
expansion_union_results_train <- read.csv(paste0("../../data/ML_out/two_class/optm/expansion_union/metrics/all_metrics_train_res.csv"))
dd_expansion_train <- read.csv(paste0("../../data/ML_out/two_class/optm/data_driven_expansion/metrics/all_metrics_train_res.csv"))
disgenet_expansion_train <- read.csv(paste0("../../data/ML_out/two_class/optm/disgenet_expansion/metrics/all_metrics_train_res.csv"))

dea_results_train$ml_input <- rep("dea", nrow(dea_results_train))
# minmax_100_results_train$ml_input <- rep("minmax 100", nrow(minmax_100_results_train))
data_driven_results_train$ml_input <- rep("data-driven", nrow(data_driven_results_train))
# minmax_30_results_train$ml_input <- rep("minmax 30", nrow(minmax_30_results_train))
# minmax_76_results_train$ml_input <- rep("minmax 76", nrow(minmax_76_results_train))
# minmax_163_results_train$ml_input <- rep("minmax 163", nrow(minmax_163_results_train))
dgn_results_train$ml_input <- rep("COPD-related", nrow(dgn_results_train))
expansion_intersection_results_train$ml_input <- rep("expansion intersection", nrow(expansion_intersection_results_train))
expansion_union_results_train$ml_input <- rep("expansion union", nrow(expansion_union_results_train))
dd_expansion_train$ml_input <- rep("data-driven expansion", nrow(dd_expansion_train))
disgenet_expansion_train$ml_input <- rep("COPD-related expansion", nrow(disgenet_expansion_train))

## TEST
dea_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/dea/metrics/all_metrics_test.csv"))
# minmax_100_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_100/metrics/all_metrics_test.csv"))
data_driven_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/data_driven/metrics/all_metrics_test.csv"))
# minmax_30_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_30/metrics/all_metrics_test.csv"))
# minmax_76_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_76/metrics/all_metrics_test.csv"))
# minmax_163_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_163/metrics/all_metrics_test.csv"))
dgn_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/disgenet/metrics/all_metrics_test.csv"))
expansion_intersection_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/expansion_intersection/metrics/all_metrics_test.csv"))
expansion_union_results_test <- read.csv(paste0("../../data/ML_out/two_class/optm/expansion_union/metrics/all_metrics_test.csv"))
dd_expansion_test <- read.csv(paste0("../../data/ML_out/two_class/optm/data_driven_expansion/metrics/all_metrics_test.csv"))
disgenet_expansion_test <- read.csv(paste0("../../data/ML_out/two_class/optm/disgenet_expansion/metrics/all_metrics_test.csv"))

dea_results_test$ml_input <- rep("dea", nrow(dea_results_test))
# minmax_100_results_test$ml_input <- rep("minmax 100", nrow(minmax_100_results_test))
data_driven_results_test$ml_input <- rep("data-driven", nrow(data_driven_results_test))
# minmax_30_results_test$ml_input <- rep("minmax 30", nrow(minmax_30_results_test))
# minmax_76_results_test$ml_input <- rep("minmax 76", nrow(minmax_76_results_test))
# minmax_163_results_test$ml_input <- rep("minmax 163", nrow(minmax_163_results_test))
expansion_intersection_results_test$ml_input <- rep("expansion intersection", nrow(expansion_intersection_results_test))
expansion_union_results_test$ml_input <- rep("expansion union", nrow(expansion_union_results_test))
dgn_results_test$ml_input <- rep("COPD-related", nrow(dgn_results_test))
data_driven_results_test$ml_input <- rep("data-driven", nrow(data_driven_results_test))
dd_expansion_test$ml_input <- rep("data-driven expansion", nrow(dd_expansion_test))
disgenet_expansion_test$ml_input <- rep("COPD-related expansion", nrow(disgenet_expansion_test))


## ----------------------------------------------------------------------------------------------------------------------------------------
# Barplot Cross-validation

## Data to plot
# results_train_76 <- rbind(dea_results_train, minmax_76_results_train)
# results_train_30 <- rbind(dgn_results_train, minmax_30_results_train)
# results_train_dd <- rbind(dea_results_train,minmax_100_results_train,data_driven_results_train)
results_train<- rbind(dgn_results_train, disgenet_expansion_train,
                              data_driven_results_train, dd_expansion_train,
                              expansion_intersection_results_train, expansion_union_results_train)
# Save the data frame with all the results by ml models and inputs for radar charts
write.csv(results_train, file = "../../data/ML_out/two_class/optm/tables/results_train.csv", row.names = FALSE)

results_train_mrmr <- rbind(minmax_30_results_train, minmax_76_results_train,
                            minmax_100_results_train, minmax_163_results_train)

data_to_plot_train <- results_train_76 %>% filter(metric =="accuracy")

## Plot
ggplot(data_to_plot_train) +
  aes(
    x = classifier,
    fill = ml_input,
    group = ml_input,
    weight = mean
  ) +
  labs(x = "ML methods",
       y = "Accuracy")+
  # scale_x_discrete(
  #   labels = c("XGB", "SVM-rad", "SVM-poly", "RF", "GLM", "KNN")) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    name = "Input genes",
    values = c(
               # `dea` = "#63ace5",
               # 'disgenet'= "#63ace5",
               # `disgenet expansion` = "#4b86b4",
               # `minmax 76` = "#e7eff6",
               # 'data driven' = "#e7eff6",
               # 'data driven expansion' = "#adcbe3",
               # 'minmax 100' = "#4b86b4",
               # 'minmax 163' = "#63ace5",
               # 'omnipath intersection' = "#76b6c4",
               # 'omnipath union' = "#2a4d69",
               `dea` = "#00968B",
               'minmax 30' = "#63ace5",
               'COPD-related'= "#325486",
               `COPD-related expansion` = "#9AAAC3",
               `minmax 76` = "#dae885",
               'data-driven' = "#19787F",
               `data-driven expansion` = "#8EBCBF",
               'minmax 100' = "#4b86b4",
               'minmax 163' = "#63ace5",
               'expansion intersection' = "#FF9999",
               'expansion union' = "#D9D9D9"
               #8EBCBF, #FF9999, #19787F, #325486, #9AAAC3
    )
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Train Bayes optimization"
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

ggsave("../../data/ML_out/two_class/optm/plots/inputs_comparison/76_train.pdf", width = 15, height = 8)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Barplot TEST

## Data to plot
results_test_76 <- rbind(dea_results_test, minmax_76_results_test)
results_test_30 <- rbind(dgn_results_test, minmax_30_results_test)
results_test_dd <- rbind(dea_results_test,minmax_100_results_test,data_driven_results_test)
results_test <- rbind(dgn_results_test, disgenet_expansion_test,
                             data_driven_results_test, dd_expansion_test,
                             expansion_intersection_results_test, expansion_union_results_test)
# Save the data frame with all the results by ml models and inputs for radar charts
write.csv(results_test, file = "../../data/ML_out/two_class/optm/tables/results_test.csv", row.names = FALSE)

results_test_mrmr <- rbind(minmax_30_results_test, minmax_76_results_test, minmax_100_results_test, minmax_163_results_test)

## PLOT
data_to_plot_test <- results_test_30 %>% filter(metric =="accuracy")

ggplot(data_to_plot_test) +
  aes(
    x = classifier,
    fill = ml_input,
    group = ml_input,
    weight = estimate
  ) +
  labs(x = "ML methods",
       y = "Accuracy")+
  # scale_x_discrete(
  #   labels = c("XGB", "SVM-rad", "SVM-poly", "RF", "GLM", "KNN")) +
  geom_bar(position = "dodge") +
  scale_fill_manual(
    name = "Input genes",
    values = c(
                  # disgenet = "#a3b899",
    #            'disgenet expansion' = "#667b68",
               dea = "#00968B",
    #            `data driven` = "#f0f1b5",
    #            'data driven expansion' = "#dae885",
    #            `minmax 30` = "#f0f1b5",
               'minmax 76' = "#dae885",
    #            'omnipath intersection' = "#68955e",
    #            'minmax 100' = "#68955e",
    #            'minmax 163' = "#454d34",
    #            'omnipath union' = "#1f270d",
               `dea` = "#63ace5",
               'minmax 30' = "#63ace5",
               'COPD-related'= "#325486",
               `COPD-related expansion` = "#9AAAC3",
               `minmax 76` = "#e7eff6",
               'data-driven' = "#19787F",
               'data-driven expansion' = "#8EBCBF",
               'minmax 100' = "#4b86b4",
               'minmax 163' = "#63ace5",
               'expansion intersection' = "#FF9999",
               'expansion union' = "#D9D9D9"
               #8EBCBF, #FF9999, #19787F, #325486, #9AAAC3
    )
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Test Bayes optimization"
  ) +
  geom_text(
    aes(label = round(estimate, 3), y = estimate),
    hjust = -0.5,
    size = 2,
    position = position_dodge(width = 0.9)
  ) +
  coord_flip() +
  theme_minimal() +
  ylim(0,1)

ggsave("../../data/ML_out/two_class/optm/plots/inputs_comparison/30_test.pdf", width = 15, height = 8)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Barplot number of genes by feature selection method
fsm_genes <- data.frame("Classifier" = c("DEA", "DisGeNet", "MinMax", "Data Driven",
                        "Data Driven OmniPath Expansion", "DisGeNet OmniPath Expansion",
                        "Union OmniPath Expansions", "Intersection OmniPath Expansions"),
                        "Genes" = c(76, 30, 100, 163, 744, 1065, 1553, 256))
fsm_genes <- data.frame("Classifier" = c("COPD-related curated", "Data-driven",
                                         "COPD-related entire",
                                         "Data-driven OmniPath", "COPD-related curated OmniPath",
                                         "Data-driven GUILDify", "COPD-related curated GUILDify",
                                         "OmniPath union", "OmniPath intersection"),
                        "Genes" = c( 30,  163, 1132, 744, 1065, 293, 161, 1553, 256))
pdf(file="../../data/Inputs_lists/inputs_plot.pdf",width = 12.3,height = 6)
ggplot(fsm_genes) +
  aes(
    x = Classifier,
    fill = Classifier,
    weight = Genes
  ) +
  labs(x = "Input selection")+
  geom_bar(position = "dodge") +
  scale_fill_manual(
    values = c(
               # DisGeNet = "#325486",
               "COPD-related entire" = "#bea9de",
               "COPD-related curated" = "#325486",
               'COPD-related curated OmniPath' = "#9AAAC3",
               'COPD-related curated GUILDify' = "#cfe2f3",
               # DEA = "#f0f1b5",
               `Data-driven` = "#19787F",
               'Data-driven OmniPath' = "#8EBCBF",
               'Data-driven GUILDify' = "#d2e7d6",
               'OmniPath intersection' = "#FF9999",
               'OmniPath union' = "#ffbaba"
    )
  ) +
  scale_x_discrete(labels = c(
    "COPD-related\ncurated",
    "COPD-related\ncurated GUILDify",
    "COPD-related\ncurated OmniPath",
    "COPD-related\nentire",
    "data-driven",
    "data-driven\nGUILDify",
    "data-driven\nOmniPath",
    "OmniPath\nintersection",
    "OmniPath\nunion"
  )) +
  theme_minimal() +
  labs(
    title = "Feature Selection Methods Number of Genes"
  )
dev.off()
## ----------------------------------------------------------------------------------------------------------------------------------------
# RadarChart dirty

# data_to_radarchat <- data_driven_results_test %>% dplyr::select(classifier, estimate,metric)
# data_to_radarchat <- filter(data_to_radarchat, metric %in% c("accuracy", "sens", "spec"))
#
# data_to_radarchat_wide <- pivot_wider(data_to_radarchat, names_from = "classifier", values_from = "estimate")
# data_to_radarchat <- as.data.frame(data_to_radarchat_wide)
# rownames(data_to_radarchat) <- data_to_radarchat_wide$metric
# colnames(data_to_radarchat) <- c("RF", "SVM-rad", "SVM-poly", "GLM", "kNN", "XGB")
# write.csv(data_to_radarchat, "../../reports/tables/")
# data_to_radarchat <- data_to_radarchat[, -1]
#
# row1 <- data.frame(rf = 1, svm_r = 1, svm_p = 1, pen_reg = 1, knn = 1, xgb = 1)
# # row2 <- data.frame(rf = 0.35, svm_r = 0.35, svm_p = 0.35, pen_reg = 0.35, knn = 0.35, xgb = 0.35)
# # row2 <- data.frame(rf = 0.65, svm_r = 0.65, svm_p = 0.65, pen_reg = 0.65, knn = 0.65, xgb = 0.65)
#
# # Unir los dataframes por columnas
# data_to_radarchat <- rbind('1'=row1,'2'=row2, data_to_radarchat)
# library(RColorBrewer)
# coul <- brewer.pal(6, "Set2")
# colors_border <- coul
#
# radarchart(data_to_radarchat,
#            #custom polygon
#            pcol=colors_border, seg = 15,
#            #custom the grid
#            cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, caxislabels=seq(0.7,1,0.1))
# # Add a legend
# legend(x=1.5, y=1.5, legend = rownames(data_to_radarchat[-c(1,2),]), bty = "n",
#        pch=20,cex=1, pt.cex=1, col = colors_border)
