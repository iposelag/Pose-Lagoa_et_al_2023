
library(easyPubMed)

search_topic1 <- "((LOC401847) AND ((Chronic Obstructive Lung Disease[MeSH Terms]) OR (Chronic Obstructive Pulmonary Diseases[MeSH Terms]) OR (COAD[MeSH Terms]) OR (COPD[MeSH Terms]) OR (Chronic Obstructive Airway Disease[MeSH Terms]) OR (Chronic Obstructive Pulmonary Disease[MeSH Terms]) OR (Airflow Obstruction, Chronic[MeSH Terms]) OR (Airflow Obstructions, Chronic[MeSH Terms]) OR (Chronic Airflow Obstructions[MeSH Terms]) OR (Chronic Airflow Obstruction[MeSH Terms])))"

data_driven <- scan("../../data/OmniPath/data_driven.txt", what = "character", sep = ",")
disgenet <- scan("../../data/raw/copd_genes/COPDDisgenet.txt", what = "character", sep = ",")
eset_genes <- scan("../../data/OmniPath/eset_genes.txt",what = "character", sep = ",")
disgenet <- intersect(disgenet, eset_genes)
input <- union(data_driven, disgenet)
dea <- scan("../../data/DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt", what = "character", sep = ",")
minmax_100 <- scan("../../data/OmniPath/minmax_100.txt", what = "character", sep = ",")

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

write.csv(summary, "../../data/OmniPath/pubmed_publications.csv")

summary <- read.csv("../../data/OmniPath/pubmed_publications.csv")

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

data_to_plot <- rbind(pubmed_minmax,pubmed_minmax_dea, pubmed_dea, pubmed_dea_dgn, pubmed_disgenet)
library(forcats)

# set the order of levels of GeneSymbol
data_to_plot$GeneSymbol <- fct_inorder(data_to_plot$GeneSymbol)
data_to_plot$input <- fct_inorder(data_to_plot$input)


pdf(file="../../reports/pubmed.pdf",width = 12,height = 10)
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
    values = c('mRMR' = "#66b2b2",
        'mRMR-DEA' = "#5f9ea0",
        `DEA` = "#008080",
        'DEA-COPD-related' = "#0086ad",
        'COPD-related'= "#325486")
  ) +
  theme(axis.text.y = element_text(size = 7))
dev.off()
