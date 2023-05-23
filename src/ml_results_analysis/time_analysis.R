#!/bin/Rscript
###############################################################################
############################ Time analysis  ################################
###############################################################################
## Comparison between time of tuning and training by input gene sets and ML model

## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(dplyr); library(ggplot2);

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data
dgn_results <- get(load("../../data/ML_out/two_class/optm/disgenet/results.Rda"))
dea_results <- get(load("../../data/ML_out/two_class/optm/dea/results.Rda"))
minmax_100_results <- get(load("../../data/ML_out/two_class/optm/mrmr/mrmr_100/results.Rda"))
data_driven_results <- get(load("../../data/ML_out/two_class/optm/data_driven/results.Rda"))
# minmax_30_results <- load(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_30/results.Rda")
# minmax_76_results <- load(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_76/results.Rda")
# minmax_163_results <- load(paste0("../../data/ML_out/two_class/optm/mrmr/mrmr_163/results.Rda")
expansion_intersection_results <- get(load("../../data/ML_out/two_class/optm/expansion_intersection/results.Rda"))
expansion_union_results <- get(load("../../data/ML_out/two_class/optm/expansion_union/results.Rda"))
dd_expansion <- get(load("../../data/ML_out/two_class/optm/data_driven_expansion/results.Rda"))
dgn_expansion <- get(load("../../data/ML_out/two_class/optm/disgenet_expansion/results.Rda"))

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot

# TIme results
dgn_time <- as.data.frame(dgn_results$times)
dgn_time$dim <- 30
dea_time <- as.data.frame(dea_results$times)
dea_time$dim <- 76
minmax_100_time <- as.data.frame(minmax_100_results$times)
minmax_100_results$dim <- 100
data_driven_time <- as.data.frame(data_driven_results$times)
data_driven_time$dim <- 163
op_intersection_time <- as.data.frame(expansion_intersection_results$times)
op_intersection_time$dim <- 256
op_union_time <- as.data.frame(expansion_union_results$times)
op_union_time$dim <- 1553
dd_expansion_time <- as.data.frame(dd_expansion$times)
dd_expansion_time$dim <- 744
dgn_expansion_time <- as.data.frame(dgn_expansion$times)
dgn_expansion_time$dim <- 1065

# ALl time results
all_results_time <- rbind(dgn_time,dea_time,minmax_100_time,data_driven_time,
                          op_intersection_time,op_union_time,dd_expansion_time,
                          dgn_expansion_time)
# Subset only the time in minutos
all_results_time <- all_results_time %>% mutate_all(~ as.numeric(gsub(" mins", "", .)))

# We need the df in long format (ggplot)
all_results_time_long <- reshape2::melt(all_results_time, id.vars = "dim")

# Plot
ggplot(all_results_time_long, aes(x = dim, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  labs(x = "Number of genes", y = "Time in minuts") +
  scale_colour_manual(name = "Clasificador",
    values = c( rf = "#6b3e26",
                svm_r = "#ffc5d9",
                svm_p = "#c2f2d0",
                knn = "#ffcb85",
                pen_reg = "#fdf5c9",
                xgb = "#ff6f69"))
ggsave("../../data/ML_out/two_class/optm/plots/time_analysis.pdf", width = 15, height = 8)
