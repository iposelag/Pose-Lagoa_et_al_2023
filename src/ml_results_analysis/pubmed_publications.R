#!/bin/Rscript
###############################################################################
###################### PubMed publications analysis  ##########################
###############################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(easyPubMed); library(forcats); library(ggplot2); library(ggbreak)

## ----------------------------------------------------------------------------------------------------------------------------------------
# PubMed search (it is only needed to run it once)

## Query to search
search_topic1 <- "((LOC401847) AND ((Chronic Obstructive Lung Disease[MeSH Terms]) OR (Chronic Obstructive Pulmonary Diseases[MeSH Terms]) OR (COAD[MeSH Terms]) OR (COPD[MeSH Terms]) OR (Chronic Obstructive Airway Disease[MeSH Terms]) OR (Chronic Obstructive Pulmonary Disease[MeSH Terms]) OR (Airflow Obstruction, Chronic[MeSH Terms]) OR (Airflow Obstructions, Chronic[MeSH Terms]) OR (Chronic Airflow Obstructions[MeSH Terms]) OR (Chronic Airflow Obstruction[MeSH Terms])))"

## Load genes sets
data_driven <- scan("../../data/Inputs_lists/data_driven.txt", what = "character", sep = ",")
disgenet <- scan("../../data/Inputs_lists/disgenet_curated.txt", what = "character", sep = ",")
eset_genes <- scan("../../data/Inputs_lists/eset_genes.txt",what = "character", sep = ",")
disgenet <- intersect(disgenet, eset_genes)
input <- union(data_driven, disgenet)
dea <- scan("../../data/DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt", what = "character", sep = ",")
minmax_100 <- scan("../../data/mrmr/random/selection/minmax_100.txt", what = "character", sep = ",")

## Run the query
summary <- data.frame(GeneSymbol = character(), "Count" = numeric(), "PubMedID" = character() )

for(i in 1:length(input)){

  search_topic <- gsub("LOC401847",input[i],search_topic1)
  my_query <- get_pubmed_ids(search_topic)
  count <- as.numeric(my_query$Count)
  if(count>0){
    t <- paste(as.character(unique(unlist(my_query$IdList))), collapse = ",")
    new <- data.frame(GeneSymbol = input[i], Count = count, PubMedID = t)
    summary <- rbind(summary, new)
  }
}

## Save pubmed publications seach
write.csv(summary, "../../data/OmniPath/pubmed_publications.csv")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Plot the number of publications

## Load pubmed publications search
summary <- read.csv("../../data/OmniPath/pubmed_publications.csv")

## Obtain the list of genes without repetitions and in a determined order
## mRMR->mRMRvsDEA->DEA->DEAvsCOPD_related->COPD_related
minmax_dea <- intersect(minmax_100, dea)
minmax_100 <- setdiff(minmax_100, minmax_dea)
dea <- setdiff(dea, minmax_dea)

dea_dgn <- intersect(dea, disgenet)
dea <- setdiff(dea, dea_dgn)
disgenet <- setdiff(disgenet, dea_dgn)

pubmed_minmax <- subset(summary, GeneSymbol %in% minmax_100)
pubmed_minmax$input <- rep("mRMR", nrow(pubmed_minmax))
pubmed_minmax$x <- c(1:nrow(pubmed_minmax))
pubmed_minmax_dea <- subset(summary, GeneSymbol %in% minmax_dea)
pubmed_minmax_dea$input <- rep("mRMR-DEA", nrow(pubmed_minmax_dea))
pubmed_minmax_dea$x <- c(1:nrow(pubmed_minmax_dea))
pubmed_dea <- subset(summary, GeneSymbol %in% dea)
pubmed_dea$input <- rep("DEA", nrow(pubmed_dea))
pubmed_dea$x <- c(1:nrow(pubmed_dea))
pubmed_dea_dgn <- subset(summary, GeneSymbol %in% dea_dgn)
pubmed_dea_dgn$input <- rep("DEA-COPD-related", nrow(pubmed_dea_dgn))
pubmed_dea_dgn$x <- c(1:nrow(pubmed_dea_dgn))
pubmed_disgenet <- subset(summary, GeneSymbol %in% disgenet)
pubmed_disgenet$input <- rep("COPD-related", nrow(pubmed_disgenet))
pubmed_disgenet$x <- c(1:nrow(pubmed_disgenet))

## Data to plot
data_to_plot <- rbind(pubmed_minmax,pubmed_minmax_dea, pubmed_dea, pubmed_dea_dgn, pubmed_disgenet)

# set the order of levels of GeneSymbol
data_to_plot$GeneSymbol <- fct_inorder(data_to_plot$GeneSymbol)
data_to_plot$input <- fct_inorder(data_to_plot$input)

## PLOT
ggplot(data_to_plot) +
  theme_minimal()+
  aes(
    x = GeneSymbol,
    fill = input,
    # group = ml_input,
    weight = as.numeric(Count)
  )+
  scale_y_break(c(150, 350),space=.2)+
  scale_y_break(c(400, 900))+
  geom_bar(position = "dodge")  +
  coord_flip()+
  scale_fill_manual(
    name = "Input genes",
    values = c('mRMR' = "#B6D7A8",
        'mRMR-DEA' = "#19787F",
        `DEA` = "#FFF1AC",
        'DEA-COPD-related' = "#0086ad",
        'COPD-related'= "#325486")
  ) +
  theme(axis.text.y = element_text(size = 7))
ggsave("../../data/Inputs_lists/pubmed.png", width = 12,height = 10)

pubmed <- read.csv("../../data/pubmed_publications.csv")

pubmed_dea <- pubmed %>% filter(GeneSymbol %in% dea) # 28/76 36,84%
sum(pubmed_dea$Count)/76
sum(pubmed_dea$Count)/28
pubmed_data_driven <- pubmed %>% filter(GeneSymbol %in% data_driven) # 45/163 27,61%
sum(pubmed_data_driven$Count)/163
sum(pubmed_data_driven$Count)/45
pubmed_disgenet <- pubmed %>% filter(GeneSymbol %in% disgenet) # 29/30
sum(pubmed_disgenet$Count)/30
sum(pubmed_disgenet$Count)/29
pubmed_mrmr <- pubmed %>% filter(GeneSymbol %in% minmax_100) # 21/100  21%
sum(pubmed_mrmr$Count)/100
sum(pubmed_mrmr$Count)/21

sum(pubmed$Count)/2499

# Cargar el paquete necesario
library(stats)

# Crear un data frame con los datos
data <- data.frame(GeneGroup = c(rep("DEA", length(dea)),
                                 # rep("Data-driven", length(data_driven)),
                                 rep("MRMR", length(mrmr)),
                                 rep("Disgenet", length(disgenet))),
                   Publications = c(dea, mrmr, disgenet))

# Normalidad
qqnorm(data[data$GeneGroup == "DEA","Publications"], main = "DEA")
qqline(data[data$GeneGroup == "DEA","Publications"])
qqnorm(data[data$GeneGroup == "MRMR","Publications"], main = "MRMR")
qqline(data[data$GeneGroup == "MRMR","Publications"])
qqnorm(data[data$GeneGroup == "Disgenet","Publications"], main = "Disgenet")
qqline(data[data$GeneGroup == "Disgenet","Publications"])

require(nortest)
by(data = data,INDICES = data$GeneGroup,FUN = function(x){ shapiro.test(x$Publications)})
# Realizar el análisis de varianza
result <- aov(Publications ~ GeneGroup, data = data)

aggregate(Publications ~ GeneGroup, data = data, FUN = mean)
aggregate(Publications ~ GeneGroup, data = data, FUN = sd)
# modelo no equilibrado
fligner.test(Publications ~ GeneGroup,data)
require(car)
leveneTest(Publications ~ GeneGroup,data,center = "median")

# Realizar el análisis de varianza
result <- aov(Publications ~ GeneGroup, data = data)
plot(result)
pairwise.t.test(x = data$Publications, g = data$GeneGroup, p.adjust.method = "holm",
                pool.sd = TRUE, paired = FALSE, alternative = "two.sided")
TukeyHSD(result)
# Imprimir los resultados
summary(result)

