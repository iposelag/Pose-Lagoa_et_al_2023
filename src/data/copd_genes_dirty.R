### Genes related with COPD ###
## ----------------------------------------------------------------------------------------------------------------------------------------
# Set seed
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Packages
library(OmnipathR);library(dplyr); library(tidyr); library(reshape2); library(ggplot2);
library(clusterProfiler); library(org.Hs.eg.db); library(pheatmap); library(sva)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load eset object
load("../../data/ExpressionSetObjects/eset.Rda")

## ----------------------------------------------------------------------------------------------------------------------------------------------
# BATCH effect correction (platform): ComBAT
batch = pData(eset)$platform_id
modcombat = model.matrix(~ dis_condition, data=pData(eset))
combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)

# Phenotype Data ExpressionObject
pData_eset <- data.frame( dis_condition = as.factor(pData(eset)$dis_condition))
pData_eset$sex = as.factor(pData(eset)$sex)
pData_eset$age = as.numeric(pData(eset)$age)
pData_eset$smoker = as.factor(pData(eset)$smoker)
pData_eset$GOLD_stage = as.factor(pData(eset)$GOLD_stage)
pData_eset$pneumocystis_colonization = as.factor(pData(eset)$pneumocystis_colonization)
pData_eset$pred_dlco = as.numeric(pData(eset)$pred_dlco)
pData_eset$emphysema = as.numeric(pData(eset)$emphysema)
pData_eset$pred_fev_prebd = as.numeric(pData(eset)$pred_fev_prebd)
pData_eset$pred_fev_postbd = as.numeric(pData(eset)$pred_fev_postbd)
pData_eset$pred_fvc_prebd = as.numeric(pData(eset)$pred_fvc_prebd)
pData_eset$pred_fvc_postbd = as.numeric(pData(eset)$pred_fvc_postbd)
pData_eset$sample_id = colnames(combat_edata)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load disgenet gene-disease associations
disgenet <- read.csv("../../data/raw/copd_genes/curated_gene_disease_associations.tsv", sep ="\t")
disgenet_copd <- filter(disgenet, diseaseName == "Chronic Obstructive Airway Disease")
dim(disgenet_copd)

# Interaction of disgenet_copd genes
all_interactions <- import_all_interactions(
  organism = 9606
) %>% 
  filter(source_genesymbol %in% disgenet_copd$geneSymbol) %>% 
  filter(curation_effort > 1)

union <- unique(c(all_interactions$target_genesymbol, disgenet_copd$geneSymbol)) 

# From the copd_genes_interactions set how many are in our eset object
dim(eset)
copd_genes <- intersect(rownames(eset), union)
length(copd_genes)
# save(copd_genes, "../../data/raw/copd_genes/copd_genes.Rda")
write(copd_genes,file="../../data/raw/copd_genes/copd_genes.txt", ncolumns = 1)
# write.table(disgenet_copd$geneSymbol, "../../data/raw/copd_genes/COPDDisgenet.txt",quote=F,sep="\t",row.names=F,col.names=F)
eset_copd_genes <- eset[copd_genes,]
dim(eset_copd_genes) # 587x301


## ----------------------------------------------------------------------------------------------------------------------------------------
# Analyze expression values of these genes

exprs_dgn <- combat_edata[row.names(combat_edata) %in% disgenet_copd$geneSymbol,]
annotation_col <- data.frame(disease_condition = pData_eset$dis_condition,
                             emphysema = pData_eset$emphysema,
                             GOLDS_stage = pData_eset$GOLD_stage,
                             sex = pData_eset$sex,
                             pred_dclo = pData_eset$pred_dlco,
                             pred_fev_prebd = pData_eset$pred_fev_prebd)
rownames(annotation_col) <- colnames(exprs(eset))
col_order <- order(pData_eset$dis_condition)

# Reorder the column annotations and the expression data
annotation_col <- annotation_col[col_order, ]
exprs_dgn <- exprs_dgn[, col_order]
pdf("../../data/ML_out/two_class/expr_copd_genes_hist.pdf")
hist(exprs_dgn)
dev.off()
pheatmap(exprs_dgn,
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(200),
         annotation_col = annotation_col, cluster_cols = FALSE, cluster_rows = TRUE,
         filename = "../../data/ML_out/two_class/expr_copd_genes.pdf",
         main = "Heatmap de expresión", fontsize_row = 6, fontsize_col = 4, fontsize = 5)
dev.off()

# gene_expr <- as.data.frame(melt(exprs(eset_disgenet)))
# 
# # Graficar el boxplot
# 
# ggplot(gene_expr, aes(x=Var1, y=value)) +
#   geom_boxplot() +
#   geom_jitter(aes(color=Var2), alpha=0.5, show.legend = FALSE) +  # capa adicional con los valores de la sample
#   # geom_line(aes(group=Var1, color=Var2), alpha=0.5, show.legend = FALSE) +  # capa adicional con las líneas que unen los puntos de los boxplots
#   xlab("Genes") +
#   ylab("Expresión") +
#   ggtitle("Boxplot de expresión génica") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----------------------------------------------------------------------------------------------------------------------------------------
# Intersections among gene sets

## DEA
# Intersection with deg (COPD vs CTRL) pvalue = 0.01 and lfc=log2(1.5)
DEA_lfc <- scan("../../data/DEA/platform_sex/CTRL_COPD_deg_symbols_0.01_0.584962500721156.txt",
                what = "character", sep = ",")
intersect(DEA_lfc, disgenet_copd$geneSymbol) # MMP1
intersect(DEA_lfc, union) 

# DEA_all <- read.csv("../../data/DEA/platform_sex/deg_copd_ctrl_COPD - CTRL.tsv", sep ="\t")
# intersect(DEA_all$gene_symbol, disgenet_copd$geneSymbol) 
# DEA_all_ordered <- DEA_all[order(DEA_all$adj.P.Val, decreasing=FALSE),]
# proba <- filter(DEA_all_ordered, DEA_all_ordered$gene_symbol %in% disgenet_copd$geneSymbol) 
exprs_dea <- combat_edata[row.names(combat_edata) %in% DEA_lfc,]
pdf("../../data/ML_out/two_class/expr_dea_hist.pdf")
hist(exprs_dea)
dev.off()

## MRMR
minmax_100 <- read.table(file = "../../data/ML_out/two_class/mRMR/mrmr/100_output.txt",
           sep = '\t', header = TRUE)
minmax_100 <- gsub(" ", "", minmax_100$Name)
# DEA_minmax <- filter(DEA_all_ordered, DEA_all_ordered$gene_symbol %in% minmax_500)

intersect(minmax_100, disgenet_copd$geneSymbol)
intersect(minmax_100, DEA_lfc)



exprs_minmax <- combat_edata[row.names(combat_edata) %in% minmax_500,]
pdf("../../data/ML_out/two_class/expr_mrmr_hist.pdf")
hist(exprs_minmax)
dev.off()

# Crear una columna con valores de edad en annotation_col
# annotation_col$emphysema <- pheno_data$emphysema
exprs_minmax <- exprs_minmax[, col_order]
pheatmap(exprs_minmax[1:70,],
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(200),
         annotation_col = annotation_col, cluster_cols = FALSE, cluster_rows = TRUE,
         filename = "../../data/ML_out/two_class/expr_mrmr_70.pdf",
         main = "Heatmap de expresión", fontsize_row = 6, fontsize_col = 4, fontsize = 5
)
dev.off()

# Crear un vector de nombres de fila vacío del mismo tamaño que el original
labels_row <- rep("", length(rownames(exprs_dea)))
list <- unique(union(disgenet_copd$geneSymbol, minmax_500))
# Asignar el nombre del gen si está en disgenet_copd$geneSymbol
labels_row <- ifelse(rownames(exprs_dea) %in% list, rownames(exprs_dea), labels_row)

## Pheatmap con genes de DEA marcando os xens de disgenet e mrmr
pheatmap(exprs_dea,
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(200),
         annotation_col = annotation_col, cluster_cols = FALSE, cluster_rows = TRUE,
         labels_row = labels_row,
         filename = "../../data/ML_out/two_class/expr_dea.pdf",
         main = "Heatmap de expresión", fontsize_row = 6, fontsize_col = 4, fontsize = 5)
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
# Enrichments analysis
# facelo con dea
enrich_result_dea <- enrichGO(gene = DEA_lfc,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.01)
# Summary results
enrich_result_dea_table <- as.data.frame(enrich_result_dea)
readr::write_tsv(
  enrich_result_dea_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_dea.tsv"
  ))

# Results dot diagram
dotplot(enrich_result_dea,
        showCategory = 50,
        title = "Análisis de enriquecimiento de genes relacionados con COPD",
        font.size = 8)
ggsave(filename = "../../data/ML_out/two_class/enrich_dea.pdf",
       width = 12, height = 10)

## DisGeNet genes enrichment analysis
enrich_result_dgn <- enrichGO(gene = disgenet_copd$geneSymbol,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          universe = rownames(eset),
                          keyType = "SYMBOL",
                          pvalueCutoff = 0.05)
# Summary results
enrich_result_dgn_table <- as.data.frame(enrich_result_dgn)
readr::write_tsv(
  enrich_result_dgn_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_dgn.tsv"
  ))

# Results dot diagram
dotplot(enrich_result_dgn,
        showCategory = 32,
        title = "Análisis de enriquecimiento de genes relacionados con COPD",
        font.size = 8)
ggsave(filename = "../../data/ML_out/two_class/enrich_dgn.pdf",
       width = 12, height = 10)

## copd_genes_interaction genes enrichment analysis
enrich_result_copd_genes <- enrichGO(gene = copd_genes,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = rownames(eset),
                              keyType = "SYMBOL",
                              pvalueCutoff = 0.01)

# Summary results
enrich_result_copd_genes_table <- as.data.frame(enrich_result_copd_genes)
readr::write_tsv(
  enrich_result_dgn_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_copd_genes.tsv"
  ))

# Results dot diagram

dotplot(enrich_result_copd_genes,
        showCategory = 50,
        title = "Análisis de enriquecimiento de genes relacionados con COPD",
        font.size = 8)
ggsave(filename = "../../data/ML_out/two_class/enrich_interactions.pdf",
       width = 12, height = 10)

## mRMR genes enrichment analysis
enrich_result_mrmr <- enrichGO(gene = minmax_100,
                               OrgDb = org.Hs.eg.db,
                               ont="BP",
                               universe = rownames(eset),
                               keyType = "SYMBOL",
                               pvalueCutoff = 0.05)

# Summary results
enrich_result_mrmr_table <- as.data.frame(enrich_result_mrmr)
readr::write_tsv(
  enrich_result_dgn_table,
  file.path(
    "../../data/Enrichment/feature_selection",
    "enrichment_mrmr.tsv"
  ))

# Results dot diagram
dotplot(enrich_result_mrmr,
        showCategory = 32,
        title = "Análisis de enriquecimiento de genes relacionados con COPD",
        font.size = 8)
ggsave(filename = "../../data/ML_out/two_class/enrich_mixmax.pdf",
       width = 12, height = 10)




# Load GWAS catalog of COPD
# gwas <- read.csv("../../data/raw/copd_genes/efotraits_EFO_0000341-associations-2023-03-14.csv")
# gwas_genes <- unique(unlist(strsplit(gwas$Mapped.gene, ", "))) # 419
# 
# ## Intersection with combined reactome genes
# # reactome_gwas <- intersect(reactome_genes, gwas_genes) # 27
# gwas_disgenet <- intersect(gwas_genes, disgenet_copd$geneSymbol) #5

# union <- unique(c(gwas_disgenet,reactome_disgenet,reactome_gwas)) # 47

# ## emphysema (Trait.s. as emphysema)
# emphysema <- gwas %>% filter(Trait.s. == "emphysema")
# 
# ## smoker (Reported trait contain word smoker)
# filtered_df <- df %>% filter(if_any(everything(), str_detect, "emphysema"))
# 
# smoker <- grepl("emphysema", gwas$Reported.trait)

# load("../../data/ExpressionSetObjects/eset.Rda")
# eset[disgenet_copd$geneSymbol,]

# # Reactome gene sets from gprofiler
# # write.table(disgenet_copd$geneSymbol, "../../data/raw/copd_genes/COPDDisgenet.txt",quote=F,sep="\t",row.names=F,col.names=F)
# reactome <- read.csv("../../data/raw/copd_genes/gProfiler_hsapiens_3-13-2023_5-54-16 PM__intersections.csv")
# reactome <- filter(reactome, source == "REAC")
# reactome_gene_set <- read.table("../../data/raw/copd_genes/reactome.v2023.1.Hs.symbols.gmt", 
#                                 sep="\t", header=FALSE, fill=TRUE, quote="", comment.char="")
# 
# # Filter enriched pathways by g:profiler (adj p value of 0.05)
# pathway_names <- paste("REACTOME_",toupper(reactome$term_name),sep="")
# ## replace blank spaces with underscores
# pathway_names <- gsub(" ", "_", pathway_names)
# ## replace hyphens with underscores
# pathway_names <- gsub("-", "_", pathway_names)
# 
# reactome_gene_set <- filter(reactome_gene_set, V1 %in% pathway_names)
# ## load inmune_system it does not appear on reactome file
# ## I wont use it finally as it has too many genes
# inmune_system <- read.table("../../data/raw/copd_genes/REACTOME_IMMUNE_SYSTEM.v2023.1.Hs.gmt", sep="\t", header=FALSE, fill=TRUE, quote="", comment.char="")
# ## load biosynthesis of protectins it does not appear on reactome file pathway with 4 genes
# ######## pendiente
# ## join all reactome gene sets
# reactome_gene_set <- rbind.fill(reactome_gene_set, inmune_system)
# 
# # unite selected columns into a single variable of unique genes
# reactome_genes <- unique(unlist(reactome_gene_set[,-c(1,2)])) 
# length(reactome_genes)
# reactome_disgenet <- intersect(reactome_genes, disgenet_copd$geneSymbol) # 15