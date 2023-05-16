#!/bin/Rscript
###############################################################################
######################## OBTAIN THE RAW DATA ##################################
###############################################################################

## REQUIRED LIBRARIES ##
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(require(biomaRt))

## Download and unzip de files
# GEOquery::getGEO("GSE47460")
# system("gzip -d *.gz")
##################################### not needed in Mare
setwd("data/raw/GSE47460")

## Create required directories
# if ("Annot"%in%list.files(".") == FALSE){dir.create("Annot")}
# if ("GPL14550"%in%list.files(".") == FALSE){dir.create("GLP14550")}
# if ("GLP6480"%in%list.files(".") == FALSE){dir.create("GLP6480")}

## Load info about the samples in each platform
platform <- data.frame(t(read.delim("PLATFORM.txt", sep = " ",
                                       header = FALSE)),
                       stringsAsFactors = FALSE)
names(platform) <- platform[1,]
platform <- platform[-1,]
platform[platform == ""] <- " "

## Files from the directory: GSE47460_RAW
all_files <- list.files(".")

## Move arrays of GLP6480
for(i in 1:length(all_files)){
  for(j in 1:length(platform[,2])){
    if(str_detect(all_files[i], platform[j,2]) == TRUE){
      file.rename(from = all_files[i], to = paste("GLP6480/",
                                                  as.name(all_files[i]), 
                                                  sep= "" ))
    }
  }
}

## Move arrays of GLP14550
for(i in 1:length(all_files)){
  for(j in 1:length(platform[,1])){
    if(str_detect(all_files[i], platform[j,1]) == TRUE){
      file.rename(from = all_files[i], to = paste("GLP14550/", 
                                                  as.name(all_files[i]), 
                                                  sep= "" ))
    }
  }
}

###############################################################################
################# OBTAIN ANNOTATIONS (using BioMart) ##########################
###############################################################################

## agilent_wholegenome_4x44k_v1
## In order to see all the Agilent tables accessible via biomaRt: 
## listAttributes(mart)[grep('agilent', tolower(listAttributes(mart)[,1])),]

# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('hsapiens_gene_ensembl', mart)
# 
# ## PLATFORM GLP6480
# annotLookup <- getBM(
#   mart = mart,
#   attributes = c(  # listAttributes(mart) to see more attributes
#     'agilent_wholegenome_4x44k_v1',
#     'wikigene_description',
#     'ensembl_gene_id',
#     'entrezgene_id',
#     'gene_biotype',
#     'external_gene_name'))
# 
# ### Save the annotations in a .tsv file
# write.table(
#   annotLookup,
#   paste0('Annot/Human_agilent_wholegenome_4x44k_v1_', 
#          gsub("-", "_", as.character(Sys.Date())),
#          '.tsv'),
#   sep = '\t',
#   row.names = FALSE,
#   quote = FALSE)
# 
# ## PLATFORM 14550
# annotLookup <- getBM(
#   mart = mart,
#   attributes = c(  # listAttributes(mart) to see more attributes
#     'agilent_sureprint_g3_ge_8x60k',
#     'wikigene_description',
#     'ensembl_gene_id',
#     'entrezgene_id',
#     'gene_biotype',
#     'external_gene_name'))
# 
# ### Save the annotations in a .tsv file
# write.table(
#   annotLookup,
#   paste0('Annot/Human_agilent_sureprint_g3_ge_8x60k_', 
#          gsub("-", "_", as.character(Sys.Date())),
#          '.tsv'),
#   sep = '\t',
#   row.names = FALSE,
#   quote = FALSE)

###############################################################################
################### OBTAIN ANNOTATIONS (using GEO) ############################
###############################################################################
# Annotation data ussing platform info from GEO
gpl_14550 <- getGEO("GPL14550")
save(gpl_14550,file="gpl_GPL14550.Rdata")
gpl_6480 <- getGEO("GPL6480")
save(gpl_6480,file="gpl_GPL6480.Rdata")

