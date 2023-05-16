#!/bin/Rscript
###############################################################################
############### OBTAIN ESET OBJECTS FOR EACH PLATFORM #########################
###############################################################################

## ----setup------------------------------------------------------------------------------------------------------------
# Set working directory
setwd("../../data/raw/GSE47460/")
set.seed(1234)


## ----------------------------------------------------------------------------------------------------------------------------------------------
# Reuired packages
library(dplyr); library(GEOquery); library(data.table); library(limma); library(stringr)


## ----------------------------------------------------------------------------------------------------------------------------------------------
# Set the platform of interest from command line
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  plt = commandArgs(trailingOnly=TRUE)
}


## ----------------------------------------------------------------------------------------------------------------------------------------------
# Create required directories
if ("ExpressionSetObjects"%in%list.files("../../") == FALSE){dir.create("../../ExpressionSetObjects")}
if (plt%in%list.files("../../ExpressionSetObjects") == FALSE){dir.create(paste("../../ExpressionSetObjects/",plt,sep=""))}

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Series Matrix from GEO
geo <- getGEO("GSE47460")
if(plt == "6480"){
  geo <- geo[[2]]
  GSM1150360 <- geo$geo_accession == "GSM1150360" # Delete sample wrongly named
  geo <- geo[,!GSM1150360]
}else{
  geo <- geo[[1]]
  GSM1150199 <- geo$geo_accession == "GSM1150199" # Delete sample wrongly named
  geo <- geo[,!GSM1150199]
}

## ---- warning=FALSE----------------------------------------------------------------------------------------------------------------------------
# Load and read raw data
GLPraw_COPD = read.maimages(dir(paste("GLP",plt, "/", sep=""),"COPD.txt"), "agilent",
                                 path = paste("GLP",plt,sep = ""), green.only = TRUE,
                                 other.columns="gIsWellAboveBG")

GLPraw_CTRL = read.maimages(dir(paste("GLP",plt, "/", sep=""), "CTRL.txt"), "agilent",
                                 path = paste("GLP",plt,sep = ""), green.only = TRUE,
                                 other.columns="gIsWellAboveBG")

# GLPraw_ILD = read.maimages(dir(paste("GLP",plt, "/", sep=""), "ILD.txt"), "agilent",
#                                 path = paste("GLP",plt,sep = ""), green.only = TRUE,
#                                 other.columns="gIsWellAboveBG")

## Merge the 3 ElistRaw objects
# GLPraw <- cbind(GLPraw_COPD, GLPraw_CTRL, GLPraw_ILD)
GLPraw <- cbind(GLPraw_COPD, GLPraw_CTRL)

# Save the data
save(GLPraw,file = paste("../../ExpressionSetObjects/",plt,"/GLPraw.Rda", sep=""))
# load(paste("../../ExpressionSetObjects/",plt,"/GLPraw.Rda", sep=""))

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Load annotation according to platform
if(plt == "6480"){
  load("gpl_GPL6480.Rdata")
  gpl <- gpl_6480
}else{
  load("gpl_GPL14550.Rdata")
  gpl <- gpl_14550
}

gpl <- gpl@dataTable@table
probe_symbol <- gpl[, c("ID", "GENE_SYMBOL")]
probe_symbol <- probe_symbol %>% 
  rename(
    ID_PROBE = ID
  )

dim(probe_symbol)
probe_symbol <- probe_symbol[match(GLPraw$genes$ProbeName, probe_symbol$ID_PROBE),]
dim(probe_symbol)
table(GLPraw$genes$ProbeName == probe_symbol$ID_PROBE) # check that annots are aligned
GLPraw$genes$ID_PROBE <- probe_symbol$ID_PROBE
GLPraw$genes$GENE_SYMBOL <- probe_symbol$GENE_SYMBOL

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Background correction
GLP_bgcorrect <- backgroundCorrect(GLPraw, method = 'normexp')

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Normalization (quantile)
GLP_bgcorrect_norm <- normalizeBetweenArrays(GLP_bgcorrect, method = "quantile")

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Check data dimensions
dim(GLP_bgcorrect_norm)

# Save the data
save(GLP_bgcorrect_norm,file = paste("../../ExpressionSetObjects/",plt,"/GLP_bgcorrect_norm.Rda", sep=""))
# load(paste("../../ExpressionSetObjects/",plt,"/GLP_bgcorrect_norm.Rda", sep=""))

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Filtering of genes: controls & noSymbol
## Controls
Control <- GLP_bgcorrect_norm$genes$ControlType==1L
## NoSymbol
NoSymbol <- GLP_bgcorrect_norm$genes$GENE_SYMBOL == "" 
GLP_bgcorrect_norm_filt <- GLP_bgcorrect_norm[!Control &
                                                !NoSymbol, ]
# Check data dimensions
dim(GLP_bgcorrect_norm_filt)

## ----------------------------------------------------------------------------------------------------------------------------------------------
## Remove annotation columns we no longer need:

GLP_bgcorrect_norm_filt$genes <- GLP_bgcorrect_norm_filt$genes[,c("ControlType", "GENE_SYMBOL")]
head(GLP_bgcorrect_norm_filt$genes)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Phenotypic vbles
phenoData <- pData(geo)
# head(phenoData)

# Rename features of interest in order to make the data manipulable
colnames(phenoData)[43] <- "emphysema"
colnames(phenoData)[44] <- "pred_dlco"
colnames(phenoData)[45] <- "pred_fev_postbd"
colnames(phenoData)[46] <- "pred_fev_prebd"
colnames(phenoData)[47] <- "pred_fvc_postbd"
colnames(phenoData)[48] <- "pred_fvc_prebd"
colnames(phenoData)[49] <- "age"
colnames(phenoData)[50] <- "dis_condition"
colnames(phenoData)[51] <- "GOLD_stage"
colnames(phenoData)[52] <- "ILD_subtype"
colnames(phenoData)[53] <- "pneumocystis_colonization"
colnames(phenoData)[54] <- "sex"
colnames(phenoData)[55] <- "smoker"

# Delete ILD from phenotypic data
phenoData$dis_condition <- as.factor(phenoData$dis_condition)
phenoData <- filter(phenoData, dis_condition == "Chronic Obstructive Lung Disease" |
                      dis_condition == "Control")
phenoData$dis_condition <- droplevels(phenoData$dis_condition)
levels(phenoData$dis_condition) <- c("COPD", "CTRL")

# ----------------------------------------------------------------------------------------------------------------------------------------------
# ExpressionSet
#experimentdata <- experimentData(geo)
metadata <- data.frame(labelDescription = colnames(phenoData))
phenotypedata = new("AnnotatedDataFrame", data = phenoData,
                    varMetadata = metadata)
metadata_f <- data.frame(labelDescription = colnames(GLP_bgcorrect_norm_filt$genes))
data_annot <- data.frame(GLP_bgcorrect_norm_filt$genes)
featuredata <- new("AnnotatedDataFrame", data = data_annot,
                   varMetadata = metadata_f)
exprs = GLP_bgcorrect_norm_filt$E
colnames(exprs) <- colnames(exprs) %>%
    str_extract('G.{9}') 
exprs <- exprs[,match(rownames(phenotypedata),colnames(exprs))]
eset_outliers <- ExpressionSet(assayData = assayDataNew(exprs = exprs),
                      phenoData =  phenotypedata,
                      # experimentData = experimentdata,
                      featureData = featuredata)
# Check data dimensions
dim(eset_outliers)

# Rare samples
table(pData(eset_outliers)$dis_condition, pData(eset_outliers)$GOLD_stage)
rare_samples <- filter(pData(eset_outliers), 
                       dis_condition == "COPD" & GOLD_stage == "0-At Risk")
table(rownames(rare_samples),rare_samples$dis_condition)
eset_outliers <- eset_outliers[,!(colnames(exprs(eset_outliers)) %in% rownames(rare_samples))]
dim(eset_outliers)
n_samples <- ncol(eset_outliers)
n_genes <- nrow(eset_outliers)

save(eset_outliers, file = paste("../../ExpressionSetObjects/",plt,"/eset_outliers.Rda", sep=""))
# load(paste("../../ExpressionSetObjects/",plt,"/eset_outliers.Rda", sep=""))

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Filtering of low quality samples (by Boxplot from arrayQualityMetrics)

results <- data.frame()

for( i in seq(1, 10000, by=50)){
  result <- read.table(paste("../../Visualization/outliers/outliers_",plt,"_",i,".txt", sep = ""))
  results <- rbind(results, result)
}

outliers <- names(subset(table(results), table(results)> 10000*0.25))

eset_outliers <- eset_outliers[,!(colnames(exprs(eset_outliers)) %in% outliers)]
dim(eset_outliers)
n_samples <- ncol(eset_outliers)
n_genes <- nrow(eset_outliers)

## ----------------------------------------------------------------------------------------------------------------------------------------------
## Filter lowly expressed genes
other <- GLP_bgcorrect_norm_filt$other$gIsWellAboveBG
colnames(other) <- colnames(other) %>%
  str_extract('G.{9}') 
other <- other[,match(rownames(phenotypedata),colnames(other))]
other <- other[,!colnames(other)
      %in% outliers]
## Is Exprs
cond <- min(table(eset_outliers$dis_condition))
IsExpr <- rowSums(other > 0) >= ceiling(cond/2)
eset_filt <- eset_outliers[IsExpr, ]
dim(eset_filt)

## ----------------------------------------------------------------------------------------------------------------------------------------------
## Compress probes into gens
expression_data <- as.data.frame(exprs(eset_filt))
expression_data$GENE_SYMBOL <- fData(eset_filt)$GENE_SYMBOL
expression_data_dt <- data.table(expression_data)
expression_data_dt <- expression_data_dt[,lapply(.SD,median),by=GENE_SYMBOL]
dim(expression_data_dt)
gene_symbol_ord <- expression_data_dt$GENE_SYMBOL
# Delete duplicated gene symbol from feature data
dup <- duplicated(fData(eset_filt)$GENE_SYMBOL)
fData(eset_filt) <- fData(eset_filt)[!dup,]
### check that gene_symbol ids are aligned
table(fData(eset_filt)$GENE_SYMBOL == gene_symbol_ord)

# Final ExpressionSet
# phenotype data is the same as before
phenotypedata <- eset_filt@phenoData
metadata_f <- data.frame(labelDescription = colnames(fData(eset_filt)))
data_annot <- data.frame(fData(eset_filt),
                          row.names = fData(eset_filt)$GENE_SYMBOL)
featuredata <- new("AnnotatedDataFrame", data = data_annot,
                   varMetadata = metadata_f)
exprs = as.matrix(as.data.frame(expression_data_dt)[, -1])
exprs <- exprs[,match(rownames(phenotypedata),colnames(exprs))]
eset <- ExpressionSet(assayData = assayDataNew(exprs = exprs),
                        phenoData =  phenotypedata,
                        # experimentData = experimentdata,
                        featureData = featuredata)
# Check data dimensions
dim(eset)


## ----------------------------------------------------------------------------------------------------------------------------------------------
# Save the data
save(eset,file = paste("../../ExpressionSetObjects/",plt,"/expression.Rda", sep=""))