# Set working directory
setwd(getSrcDirectory(function(){})[1]) # Not working on Rstudio
# Set seed reproducibility
set.seed(1234)

## --------------------------------------------------------------------------------------------------------------------
# Required packages
library(Biobase); library(limma); library(sva);

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Read dea cutoffs an input eset object from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Please provide one arguments.")
} else if (length(args)==3) {
  esetobject <- args[1]
  pv_c <- args[2]
  lfc_c <- args[3]
}

pv <- as.numeric(pv_c)
lfc <- as.numeric(lfc_c)
esetobject <- "eset"
if(esetobject == "eset_several_class"){
  classif <- "three_class"
}else{
  classif <- "two_class"
}

## --------------------------------------------------------------------------------------------------------------------
# Load data
eset <- get(load(paste("../../data/ExpressionSetObjects/",esetobject,".Rda",sep = "" )))
samples_train <- read.table(paste0("../../data/ML_out/",classif,"/samples_train.txt"))
# Check data dimensions
print("Dimensions eset: platforms combination")
dim(eset) # 16235x318 all_Samples
eset <- eset[,samples_train$V1]
dim(eset)
# 16235x301 filtered
n_samples <- ncol(eset)
n_genes <- nrow(eset)

## --------------------------------------------------------------------------------------------------------------------
# DEA
## --------------------------------------------------------------------------------------------------------------------

# Define the function two class
## It requires:
##  groups: levels of disease_condition
##  comb: combination of these levels
##  pv: pvalue
##  lfc: log fold change
myDEA2 <- function(eset, groups, comb, pv = 0.05, lfc = 0) {

  # Extract gene symbols
  genes <- data.frame(gene_symbol = fData(eset)$GENE_SYMBOL)

  # Create model matrix
  mod <- model.matrix(~ 0 + dis_condition + platform_id + sex, data = pData(eset))
  colnames(mod) <- c(groups, "GPL6480", "Female")

  # Fit the linear model
  fit <- lmFit(exprs(eset), mod)

  # Set up the contrasts
  contrast.matrix <- makeContrasts(
    comb,
    levels = mod
  )

  # Perform differential expression analysis
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  results <- decideTests(fit1, adjust.method = "fdr", p.value = pv, lfc = lfc)

  # Add gene symbols to the results
  fit1$genes <- genes

  # Plot a Venn diagram of the differentially expressed genes
  pdf(paste0("../../data/DEA/platform_sex/venn_", paste0(groups, collapse = "_"),"_",pv,"_",round(lfc,2), ".pdf"))
  vennDiagram(results, include = c("up", "down"))
  dev.off()

  # Print a summary of the results
  print(summary(results))

  # Top differentially expressed genes between groups:

  degenes <- c()
  for(i in 1:length(comb)){
    pv1 <- fit1$p.value[,comb[i]]        # p value
    lfc1 <- fit1$coefficients[,comb[i]]  # log fold change
    n <- sum(p.adjust(pv1, method="fdr")<pv & abs(lfc1) > lfc)
    print(n)
    top_deg <- topTable(fit1, coef=i,  number = n, adjust.method = "fdr", p.value = pv, lfc = lfc, sort.by = "p")
    # print(head(top_deg))
    # save all diff expr gene values
    deg_all <- topTable(fit1, coef=i,  number = nrow(eset), adjust.method = "fdr", sort.by = "p")
    readr::write_tsv(
      deg_all,
      file.path(
        "../../data/DEA/platform_sex",
        paste("deg_",comb[i],".tsv",sep="")
      ))
    # save results
    readr::write_tsv(
      top_deg,
      file.path(
        "../../data/DEA/platform_sex",
        paste("deg_",comb[i],"_",pv,"_",round(lfc,2),".tsv",sep="")
      ))
    # png()
    # hist_deg <- hist(pv_c, main="Distribution of expected p-values",
    #                  las=1, ylim = c(0, 8000))
    # dev.off()
    plotMD(fit1, coef=1, status=results)
    minmax100 <- scan("../../data/OmniPath/minmax_100.txt", what = "character", sep = ",")
    data_driven_genes <- rownames(fit1)[rownames(fit1) %in% minmax100]
    #Add labels to the top genes
    points(x = fit1$Amean[data_driven_genes], y = fit1$coefficients[,1][data_driven_genes], pch = 16, col = "#36D9AD")
    text(x = fit1$Amean[data_driven_genes], y = fit1$coefficients[,1][data_driven_genes], labels = data_driven_genes, col="#36D9AD",cex=0.6,pos=3, font = 2)
    ggplot2::ggsave(file.path("../../DEA/platform_sex", paste("plotmd_copd_ctrl_",pv_c,".png",sep="")),
                    plot = plotmd
    )
    degenes <- c(degenes,top_deg$gene_symbol)

  }

  degenes <- unique(degenes)
  print(degenes)
  write(degenes,file=paste0("../../data/DEA/platform_sex/",paste0(groups, collapse = "_"),"_deg_symbols_",pv,"_",lfc,".txt"), ncolumns = 1)

  # Return the results
  return(top_deg)
}

## --------------------------------------------------------------------------------------------------------------------
# Define the function (three class)

## It requires:
##  groups: levels of disease_condition
##  comb: combination of these levels
##  pv: pvalue
##  lfc: log fold change
myDEA3 <- function(eset, groups, comb, pv = 0.05, lfc = 0) {

  # Extract gene symbols
  genes <- data.frame(gene_symbol = fData(eset)$GENE_SYMBOL)

  # Create model matrix
  mod <- model.matrix(~ 0 + dis_condition + platform_id + sex, data = pData(eset))
  colnames(mod) <- c(groups, "GPL6480", "Female")

  # Fit the linear model
  fit <- lmFit(exprs(eset), mod)

  # Set up the contrasts
  contrast.matrix <- makeContrasts(
    COPD_Severe-CTRL_At_Risk,
    COPD_Severe-COPD_No_Severe,
    COPD_No_Severe-CTRL_At_Risk,
    levels = mod
  )

  # Perform differential expression analysis
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  results <- decideTests(fit1, adjust.method = "fdr", p.value = pv, lfc = lfc)

  # Add gene symbols to the results
  fit1$genes <- genes

  # Plot a Venn diagram of the differentially expressed genes
  pdf(paste0("../../data/DEA/platform_sex/venn_", paste0(groups, collapse = "_"),"_",pv,"_",round(lfc,2), ".pdf"))
  vennDiagram(results, include = c("up", "down"))
  dev.off()

  # Print a summary of the results
  print(summary(results))

  # Top differentially expressed genes between groups:

  degenes <- c()
  for(i in 1:length(comb)){
    pv1 <- fit1$p.value[,comb[i]]        # p value
    lfc1 <- fit1$coefficients[,comb[i]]  # log fold change
    n <- sum(p.adjust(pv1, method="fdr") < pv & abs(lfc1) > lfc)
    print(n)
    top_deg <- topTable(fit1, coef=i,  number = n, adjust.method = "fdr", p.value = pv, lfc = lfc, sort.by = "p")
    # print(head(top_deg))
    deg_all <- topTable(fit1, coef=i,  number = nrow(eset), adjust.method = "fdr", sort.by = "p")
    readr::write_tsv(
      deg_all,
      file.path(
        "../../data/DEA/platform_sex",
        paste("deg_",comb[i],".tsv",sep="")
      ))
    # save results
    readr::write_tsv(
      top_deg,
      file.path(
        "../../data/DEA/platform_sex",
        paste("deg_",comb[i],"_",pv,"_",round(lfc,2),".tsv",sep="")
      ))
    # png()
    # hist_deg <- hist(pv_c, main="Distribution of expected p-values",
    #                  las=1, ylim = c(0, 8000))
    # dev.off()
    # plotmd <- plotMD(fit1, coef=1, status=results1)
    # ggplot2::ggsave(file.path("../../DEA/platform_sex", paste("plotmd_copd_ctrl_",pv_c,".png",sep="")),
    #                 plot = plotmd
    # )
    degenes <- c(degenes,top_deg$gene_symbol)

  }

  degenes <- unique(degenes)
  print(degenes)
  write(degenes,file=paste0("../../data/DEA/platform_sex/",paste0(groups, collapse = "_"),"_deg_symbols_",pv,"_",lfc,".txt"), ncolumns = 1)

  # Return the results
  return(top_deg)
}

## --------------------------------------------------------------------------------------------------------------------
# Run DEA

if(esetobject == "eset_several_class"){
  eset$dis_condition <- relevel(factor(eset$dis_condition), ref="CTRL_At_Risk")
  groups <- levels(eset$dis_condition)
  comb <- c("COPD_Severe - CTRL_At_Risk",
            "COPD_Severe - COPD_No_Severe",
            "COPD_No_Severe - CTRL_At_Risk")

  myDEA3(eset, groups, comb, pv = 0.01, lfc = log2(1.5))
}else{
  eset$dis_condition <- relevel(factor(eset$dis_condition), ref="CTRL")
  groups <- levels(eset$dis_condition)
  comb <- c("COPD - CTRL")
  myDEA2(eset, groups, comb, pv = 0.01, lfc = log2(1.5) )
}


# Crear una function de dea para 2 grupos y para 3
# Return the resulting vector
