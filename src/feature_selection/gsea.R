#!/bin/Rscript
###############################################################################
################### Gene Set Enrichment Analysis #############################
###############################################################################
## command: Rscript gsea.R
## output: gene set expression analysis

# Set seed reproducibility
set.seed(1234)

## --------------------------------------------------------------------------------------------------------------------
# Required packages
library(Biobase); library(ReactomePA); library(ggplot2); library(dplyr); library(org.Hs.eg.db);
library(enrichplot)

## --------------------------------------------------------------------------------------------------------------------
# Required directories
if ("Enrichment"%in%list.files("../../data/") == FALSE){
  dir.create("../../data/Enrichment/", showWarnings = FALSE)}
# we will adjust by platform and sex
if ("GSEA"%in%list.files("../../data/Enrichment/") == FALSE){
  dir.create("../../data/Enrichment/GSEA", showWarnings = FALSE)}

## --------------------------------------------------------------------------------------------------------------------
# Set working directory
setwd("../../data/Enrichment/GSEA/")

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Read dea cutoffs an input eset object from command line
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 1) {
#   stop("Please provide one arguments.")
# } else if (length(args) == 1) {
#   deg_file <- args[1]
# }

a <- 'deg_COPD - CTRL'

## --------------------------------------------------------------------------------------------------------------------
# Load data

deg <- read.table("../../data/DEA/platform_sex/deg_COPD - CTRL.tsv", sep = "\t", header = TRUE)

## --------------------------------------------------------------------------------------------------------------------
# Change gene symbol to entrezID

# Example gene symbols
gene_symbols <- deg$gene_symbol

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", column = "ENTREZID")
entrez_ids <- as.data.frame(entrez_ids)
entrez_ids$gene_symbol <- rownames(entrez_ids)

# Add entrez_ids to deg dataframe
deg <- merge(deg, entrez_ids, by = "gene_symbol")

## --------------------------------------------------------------------------------------------------------------------
# GSEA
## --------------------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------------------
# GSEA: more
filtered_deg_mapped_df <- deg %>%
  # Sort so that the highest absolute values of the t-statistic are at the top
  dplyr::arrange(dplyr::desc(abs(t)))

# Let's create a named vector ranked based on the t-statistic values
t_vector <- filtered_deg_mapped_df$t
names(t_vector) <- filtered_deg_mapped_df$entrez_ids

# We need to sort the t-statistic values in descending order here
t_vector <- sort(t_vector, decreasing = TRUE)

# Run GSEA
gsea_results <- gsePathway(
  geneList = t_vector, # Ordered ranked gene list
  organism = "human",
  minGSSize = 1, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  eps = 0,
  pvalueCutoff = 0.05, # p-value cutoff
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  by = "fgsea"
  )

# We can access the results from our gseaResult object using `@result`
print(head(gsea_results@result))
gsea_result_df <- data.frame(gsea_results@result)

# Save results
readr::write_tsv(
  gsea_result_df,
  file.path("../../data/Enrichment/GSEA",
            paste("gsea_results_",a,".tsv",sep="")
  ))

# PLot 10 most enriched pathways (up&down)
#enrichMap(gsea_results)
d_plot <- dotplot(gsea_results, showCategory=10, split=".sign", font.size=6) +
  facet_grid(.~.sign)

ggplot2::ggsave(file.path("../../data/Enrichment/GSEA", paste("dot_plot_",a,".png",sep="")),
                  plot = d_plot)

# Plot emapplot
#' Enrichment map organizes enriched terms into a network with edges connecting
#' overlapping gene sets. In this way, mutually overlapping gene sets are tend
#' to cluster together, making it easy to identify functional modules.
x <- pairwise_termsim(gsea_results)
ep <- emapplot(x, cex_label_category = 0.4, max.overlaps = Inf)
ggplot2::ggsave(file.path("../../data/Enrichment/GSEA", paste("emma_plot_",a,".png",sep="")),
                  plot = ep)

# Plot cnetplot
#' The cnetplot depicts the linkages of genes and biological concepts (e.g. GO
#' terms or KEGG pathways) as a network (helpful to see which
#' genes are involved in enriched pathways and genes that may belong to
#' multiple annotation categories)
log_fc_vector <- filtered_deg_mapped_df$logFC
names(log_fc_vector) <- filtered_deg_mapped_df$entrez_ids
cp <- cnetplot(gsea_results, categorySize="pvalue", foldChange=log_fc_vector,
               cex_label_gene = 0.2, cex_label_category = 0.5, max.overlaps = Inf)
ggplot2::ggsave(file.path("../../data/Enrichment/GSEA", paste("cnet_plot_",a,".png",sep="")),
                  plot = cp)

gsea_result_df <- gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(n = 140, order_by = NES)

gene_set_id <- gsea_result_df[1,1]

# PLOT: most positive pathway according to NES plot
most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = gene_set_id,
  title = gene_set_id,
  color.line = "#0d76ff"
)

most_positive_nes_plot

# To save the plot
# ggplot2::ggsave(file.path("../../data/Enrichment/GSEA", paste("enr_pos_platform_sex_",gene_set_id,".png",sep="")),
#                   plot = most_positive_nes_plot)
