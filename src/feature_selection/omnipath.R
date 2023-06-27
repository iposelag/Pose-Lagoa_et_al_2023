#!/bin/Rscript
###############################################################################
############# COPD-related, data-driven and expansions #######################
###############################################################################
## command: Rscript omnipath.R
## output: Lists of genes and enrichment analysis

### Genes related with COPD ###
## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(OmnipathR);library(dplyr); library(tidyr); library(reshape2); library(ggplot2);
library(clusterProfiler); library(org.Hs.eg.db); library(pheatmap); library(sva)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("Inputs_lists"%in%list.files("../../data/") == FALSE){
  dir.create("../../data/Inputs_lists/", showWarnings = FALSE)}
if ("feature_selection"%in%list.files("../../data/Enrichment") == FALSE){
  dir.create("../../data/Enrichment/feature_selection", showWarnings = FALSE)}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

## Load eset object
load("../../data/ExpressionSetObjects/eset.Rda")
write(rownames(eset),file="../../data/Inputs_lists/eset_genes.txt", ncolumns = 1)
length(rownames(eset))

## Load disgenet gene-disease associations curated database
disgenet <- read.csv("../../data/raw/disgenet_tables/curated_gene_disease_associations.tsv", sep ="\t")
disgenet_copd <- filter(disgenet, diseaseName == "Chronic Obstructive Airway Disease")
dim(disgenet_copd)
write(disgenet_copd$geneSymbol,file="../../data/Inputs_lists/disgenet_curated.txt", ncolumns = 1)
write(intersect(disgenet_copd$geneSymbol,rownames(eset)),file="../../data/GUILDify/disgenet_curated_guildify.txt",  sep = ";", ncolumns = length(data_driven))

## Load disgenet gene-disease associations entire list (already filtered by coad)
disgenet_entire_list <- read.csv("../../data/raw/disgenet_tables/C0024117_disease_gda_summary.tsv", sep ="\t")
dim(disgenet_entire_list)
write(disgenet_entire_list$Gene,file="../../data/Inputs_lists/disgenet_entire_list.txt", ncolumns = 1)

## Load evidencies to check the curated list extracted from disgenet
# CURATED: GDAs from UniProt, PsyGeNET, Orphanet, the CGI, CTD (human data), ClinGen, and the Genomics England PanelApp
disgenet_evidences <- read.csv("../../data/raw/disgenet_tables/C0024117_disease_gda_evidences.tsv", sep = "\t")
disgenet_evidences_curated_new <- subset(disgenet_evidences, Original_DB == "CTD_human" | 
                                           Original_DB == "GENOMICS_ENGLAND")
intersect(unique(disgenet_evidences_curated_new$Gene), disgenet_copd$geneSymbol) ## the curated genes are the same as my original table

## Loadd DEg (COPD vs CTRL) pvalue = 0.01 and lfc=log2(1.5)
DEA_lfc <- scan("../../data/DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt",
                what = "character", sep = ",")

## Load MinMax
minmax <- scan(file = "../../data/mrmr/random/selection/minmax_100.txt",
                   what = "character", sep = ",")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Relations among dgn, deg, minmax
dea_dgn <- intersect(DEA_lfc, disgenet_copd$geneSymbol) # 1
dea_dgn_entire_list <- intersect(DEA_lfc, disgenet_entire_list$Gene) # 14
dea_minmax <- intersect(DEA_lfc, minmax) # 13
dgn_minmax <- intersect(disgenet_copd$geneSymbol, minmax) # 0
dgn_entire_list_minmax <- intersect(DEA_lfc, minmax) # 13

## ----------------------------------------------------------------------------------------------------------------------------------------
# Expansion with Guildify
## Data Driven 1% == 130 +98 = 228
data_driven_guildify <- read.table("../../data/GUILDify/data_driven_guild_scores.txt", head = TRUE, sep ="\t")
data_driven_guildify_list <- data_driven_guildify$Gene.Symbol[98:228]
data_driven_guildify_list <- union(data_driven, data_driven_guildify_list)
length(data_driven_guildify_list) #193
write(data_driven_guildify_list,file="../../data/Inputs_lists/data_driven_guildify.txt", ncolumns = 1)

## Disgenet curated 1% == 131 + 28 = 159
disgenet_curated_guildify <- read.table("../../data/GUILDify/disgenet_curated_guild_scores.txt", head = TRUE, sep ="\t")
disgenet_curated_guildify_list <- disgenet_curated_guildify$Gene.Symbol[28:159]
disgenet_curated_guildify_list <- union(intersect(disgenet_copd$geneSymbol,rownames(eset)),
                                        disgenet_curated_guildify_list)
length(disgenet_curated_guildify_list) 
write(disgenet_curated_guildify_list,file="../../data/Inputs_lists/disgenet_curated_guildify.txt", ncolumns = 1)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Expansion with OmniPath

## DGN
expansion_dgn <- import_all_interactions(
  organism = 9606,
  directed  = 'no'
) %>%
  filter(source_genesymbol %in% disgenet_copd$geneSymbol |
           target_genesymbol %in% disgenet_copd$geneSymbol ) %>%
  filter(curation_effort > 1)

## Number of times that each of the genes of my list appears
tt <- table(c(expansion_dgn$source_genesymbol,expansion_dgn$target_genesymbol))
sort(tt[disgenet_copd$geneSymbol])

## Union dgn + expanssion
expansion_dgn_genes <- union(expansion_dgn$target_genesymbol,expansion_dgn$source_genesymbol)
length(expansion_dgn_genes)
## Dgn genes that are not in expanssion
diff_dgn <- setdiff(disgenet_copd$geneSymbol, expansion_dgn_genes)
print(diff_dgn)
length(diff_dgn)
## Add the dgn genes that are not in the expanssion
expansion_dgn_genes <- union(expansion_dgn_genes, disgenet_copd$geneSymbol)
length(expansion_dgn_genes)
## Genes in my eset object
expansion_dgn_genes_eset <- intersect(rownames(eset), expansion_dgn_genes)
length(expansion_dgn_genes_eset)
write(expansion_dgn_genes_eset,file="../../data/OmniPath/expansion_dgn_genes.txt", ncolumns = 1)

## DEA + MinMax (data related)
length(dea_minmax) # 13
data_driven <- union(DEA_lfc, minmax)
length(data_driven) # 166
write(data_driven,file="../../data/Inputs_lists/data_driven.txt", ncolumns = 1)
write(data_driven,file="../../data/GUILDify/data_driven_guildify.txt", sep = ";", ncolumns = length(data_driven))

expansion_data_driven <- import_all_interactions(
  organism = 9606,
  directed  = 'no'
) %>%
  filter(source_genesymbol %in% data_driven |
           target_genesymbol %in% data_driven) %>%
  filter(curation_effort > 1)

## Number
tt <- table(c(expansion_data_driven$source_genesymbol,expansion_data_driven$target_genesymbol))
sort(tt[data_driven])

# Union dgn + expansion
expansion_data_driven_genes <- union(expansion_data_driven$target_genesymbol,expansion_data_driven$source_genesymbol)
length(expansion_data_driven_genes)
# Data related genes that are not in expanssion
diff_data_driven <- setdiff(data_driven, expansion_data_driven_genes)
print(diff_data_driven)
length(diff_data_driven)
# Add the data related genes that are not in the expanssion
expansion_data_driven_genes <- union(expansion_data_driven_genes, data_driven)
length(expansion_data_driven_genes)
## Genes in my eset object
expansion_data_driven_genes_eset <- intersect(rownames(eset), expansion_data_driven_genes)
length(expansion_data_driven_genes_eset)
write(expansion_data_driven_genes_eset,file="../../data/OmniPath/expansion_data_driven_genes.txt", ncolumns = 1)

## Intersection among both expansions
expansion_intersection <- intersect(expansion_data_driven_genes, expansion_dgn_genes)
length(expansion_intersection)
## Union among both expansions
expansion_union <- union(expansion_data_driven_genes, expansion_dgn_genes)
length(expansion_union)

eset_expansion_intersection <- intersect(rownames(eset), expansion_intersection)
length(eset_expansion_intersection)
write(eset_expansion_intersection,file="../../data/OmniPath/eset_expansion_intersection.txt", ncolumns = 1)
eset_expansion_union <- intersect(rownames(eset), expansion_union)
length(eset_expansion_union)
write(eset_expansion_union,file="../../data/OmniPath/eset_expansion_union.txt", ncolumns = 1)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Enrichment analysis

## DEA ##
enrich_result_dea <- enrichGO(gene = DEA_lfc,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.05)
## Summary results
enrich_result_dea_table <- as.data.frame(enrich_result_dea)

## Save results
readr::write_tsv(
  enrich_result_dea_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_dea.tsv"
  ))


## MinMax ##
enrich_result_minmax <- enrichGO(gene = minmax,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.05)
## Summary results
enrich_result_minmax_table <- as.data.frame(enrich_result_minmax)

## Save results
readr::write_tsv(
  enrich_result_minmax_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_minmax.tsv"
  ))

## DGN ##
enrich_result_dgn <- enrichGO(gene = disgenet_copd$geneSymbol,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.05)
## Summary results
enrich_result_dgn_table <- as.data.frame(enrich_result_dgn)

## Save results
readr::write_tsv(
  enrich_result_dgn_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_dgn.tsv"
  ))

## DGN EXPANSION ##
enrich_result_dgn_expansion <- enrichGO(gene = expansion_dgn_genes_eset,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.05)
## Summary results
enrich_result_dgn_expansion_table <- as.data.frame(enrich_result_dgn_expansion)

## Save results
readr::write_tsv(
  enrich_result_dgn_expansion_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_dgn_expansion.tsv"
  ))

## DATA DRIVEN ##
enrich_result_dd <- enrichGO(gene = data_driven,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.05)
## Summary results
enrich_result_dd_table <- as.data.frame(enrich_result_dd)

## Save results
readr::write_tsv(
  enrich_result_dd_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_data_driven.tsv"
  ))

## DATA DRIVEN EXPANSION ##
enrich_result_dd_expansion <- enrichGO(gene = expansion_data_driven_genes,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.05)
## Summary results
enrich_result_dd_expansion_table <- as.data.frame(enrich_result_dd_expansion)

## Save results
readr::write_tsv(
  enrich_result_dd_expansion_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_data_driven_expansion.tsv"
  ))

## INTERSECTION EXPANSION ##
enrich_result_intersection_expansion <- enrichGO(gene = eset_expansion_intersection,
                                       OrgDb = org.Hs.eg.db,
                                       ont = "BP",
                                       universe = rownames(eset),
                                       keyType = "SYMBOL",
                                       pvalueCutoff = 0.05)
## Summary results
enrich_result_intersection_expansion_table <- as.data.frame(enrich_result_intersection_expansion)

## Save results
readr::write_tsv(
  enrich_result_intersection_expansion_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_intersection_expansion.tsv"
  ))

## UNION EXPANSION ##
enrich_result_union_expansion <- enrichGO(gene = eset_expansion_union,
                                                 OrgDb = org.Hs.eg.db,
                                                 ont = "BP",
                                                 universe = rownames(eset),
                                                 keyType = "SYMBOL",
                                                 pvalueCutoff = 0.05)
## Summary results
enrich_result_union_expansion_table <- as.data.frame(enrich_result_union_expansion)

## Save results
readr::write_tsv(
  enrich_result_union_expansion_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_union_expansion.tsv"
  ))

## Intersections
dea_dgn_pathways <- intersect(enrich_result_dea_table$ID, enrich_result_dgn_table$ID)
print("DEA vs COPD-related pathways intersections:")
print(dea_dgn_pathways)
dea_minmax_pathways <- intersect(enrich_result_dea_table$ID, enrich_result_minmax_table$ID)
print("DEA vs mRMR pathways intersections:")
print(dea_minmax_pathways)
dgn_minmax_pathways <- intersect(enrich_result_dgn_table$ID, enrich_result_minmax_table$ID)
print("mRMR vs COPD-related pathways intersections:")
print(dgn_minmax_pathways)
