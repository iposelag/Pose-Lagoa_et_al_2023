#!/bin/Rscript
###############################################################################
################# randomize sample labels for mrmr ###########################
###############################################################################
## This file generates 1000 new dataset with shuffled sample labels for randomizations
## Then, it generates the null distribution of the scores obatained in mrmr and
## computes the treshold x / P(xin[x,1]<= 0.01

library(dplyr); library(ggplot2);
## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot top 500 mRMR scoresof my original data
# Read the file into a data frame
data <- read.table("../../../data/ML_out/two_class/mRMR/mrmr/output_filtered.txt", header = TRUE)
data_100 <- data[1:100,]
# Create a line plot of the scores
plot(data$Score, type = "l", xlab = "Order", ylab = "Score", main = "Scores Plot")
abline(h = 0.03, col = "red", lty = 2)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Create dataset altering the disease condition level of samples for randomize mRMR scores
# copd_ctrl_train <- read.csv("../../data/mrmr/copd_ctrl_train.csv", stringsAsFactors = F,sep=",")
#
# for(i in 389:1000){
#
  # # Shuffle the column indices
  # shuffled_names <- sample(copd_ctrl_train[,1])
  #
  # # Reorder the columns in the data frame
  # copd_ctrl_train[,1] <- shuffled_names
  #
  # # Write the shuffled data to a new CSV file
  # write.csv(copd_ctrl_train, paste0("../../data/mrmr/random/copd_ctrl_train_",i, ".csv"), row.names = FALSE)
#
# }

## ----------------------------------------------------------------------------------------------------------------------------------------
# Analyze randomizations

## Distribuciones de las randomizaciones

all_random <- data.frame()

# set the folder path
folder_path <- "../../../mrmr/random/output_filtered/"

# get the file names that match the pattern
file_names <- list.files(folder_path, pattern = "\\d+\\_output.txt")

# read in the files
for( i in 1:length(file_names)){
  # print(i)
  iteration <- data.frame(iteration  = rep(i,500))
  mrmr_results <- read.table(paste0(folder_path,file_names[i]), header = TRUE)

  score <- data.frame(score = mrmr_results$Score)

  random <- cbind(iteration, score)
  all_random <- rbind(random, all_random)
}

# Check if our data is normally distributed
ks.test(all_random$score, "pnorm", mean = mean(all_random$score), sd = sd(all_random$score))
## pvalue to small -> rechazamos H0 de normalidad

## Asumiendo q la distribucion no es normal
sorted_data <- sort(all_random$score, decreasing = TRUE)
n <- length(sorted_data)
index <- ceiling(0.01 * n)
a <- sorted_data[index]

## Genes in our list surpassing the treshold
mrmr_genes <- filter(data, Score > a) #103
write.csv(mrmr_genes, file="../../../mrmr/random/selection/mrmr_genes_103.csv", row.names = FALSE)
write(mrmr_genes$Name,file="../../../mrmr/random/selection/minmax_103.txt", ncolumns = 1)

## Prob
ecdf_function <- ecdf(all_random$score)
prob <- ecdf_function(1) - ecdf_function(a)
print(prob)

# Plot distribution and treshold obtained
plot <- ggplot(all_random, aes(x = score)) +
  geom_density(fill = "gray", color = "black") +
  geom_vline(xintercept = a)+
  theme_minimal()

new_plot <- plot +
  geom_text(aes(x = 0.1, y = 0.9, label = "Prob < 0.01"),
            color = "red", vjust = -0.5, size = 4)

print(new_plot)
