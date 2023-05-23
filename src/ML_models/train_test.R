#!/bin/Rscript
###############################################################################
################## GENERATION TRAIN AND TEST SETS ############################
###############################################################################
## Generate train and test sets for eset (COPD vs CTRL) and eset_several_class (CTRL vs SEV vs MSEV)
## command: Rscript train_test.R esetobject (Rscript train_test.R eset)
## output: differential expression analysis results

# Set seet (reproducibility)
set.seed(1234)

## --------------------------------------------------------------------------------------------------------------------
# Required packages
library(Biobase); library(sva); library(rsample)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Read dea cutoffs an input eset object from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide one arguments.")
} else if (length(args)==1) {
  esetobject <- args[1]
}

if(esetobject == "eset_several_class"){
  classif <- "three_class"
}else{
  classif <- "two_class"
}

eset <- get(load(paste("../../data/ExpressionSetObjects/",esetobject,".Rda",sep = "" )))

## ----------------------------------------------------------------------------------------------------------------------------------------------
# Required directories
if ("ML_out"%in%list.files("../../data/") == FALSE){
  dir.create("../../data/ML_out/", showWarnings = FALSE)}
if (classif%in%list.files("../../data/ML_out/") == FALSE){
  dir.create(paste0("../../data/ML_out/",classif), showWarnings = FALSE)}
if ("mRMR"%in%list.files(paste0("../../data/ML_out/", classif)) == FALSE){
  dir.create(paste0("../../data/ML_out/",classif, "/mRMR"), showWarnings = FALSE)}

## ----------------------------------------------------------------------------------------------------------------------------------------------
# BATCH effect correction (platform): ComBAT
batch = pData(eset)$platform_id
modcombat = model.matrix(~ dis_condition, data=pData(eset))
combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)
# plotMDS(combat_edata)

# Data frame of interest with ComBat data
combat_edata_ph <- as.data.frame(t(combat_edata))
dim(combat_edata_ph)

combat_edata_ph$dis_condition <- as.factor(pData(eset)$dis_condition)
combat_copd_ctrl <- combat_edata_ph
dim(combat_copd_ctrl)
combat_copd_ctrl = combat_copd_ctrl[ ,
                                     c("dis_condition",names(combat_copd_ctrl)[
                                       names(combat_copd_ctrl) != "dis_condition"])]
dim(combat_copd_ctrl)
combat_gold_stage <- data.frame(GOLD_stage = as.factor(pData(eset)$GOLD_stage))
combat_gold_stage = cbind(combat_gold_stage, combat_edata_ph)
combat_gold_stage = combat_gold_stage[ ,
                                     names(combat_gold_stage)[
                                       names(combat_gold_stage) != "dis_condition"]]
dim(combat_gold_stage)

# write.table(combat_copd_ctrl, "../../data/ML_out/two_class/copd_ctrl_expression.csv",quote=F,sep="\t",row.names=F,col.names=T)
# write.table(combat_gold_stage, "../../data/ML_out/two_class/gold_stage_expression.csv",quote=F,sep="\t",row.names=F,col.names=T)

if(esetobject == "eset_several_class"){
  combat_copd_ctrl$dis_condition <- relevel(combat_copd_ctrl$dis_condition,
                                            ref="COPD_Severe")
}else{
  combat_copd_ctrl$dis_condition <- relevel(combat_copd_ctrl$dis_condition,
                                            ref="COPD")
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Split on Train and Test for ML models
## ----------------------------------------------------------------------------------------------------------------------------------------

# Split data into train and test set
copd_ctrl_split <- initial_split(combat_copd_ctrl, strata = dis_condition)
copd_ctrl_split

# Train
copd_ctrl_train <- training(copd_ctrl_split)
dim(copd_ctrl_train)
write(rownames(copd_ctrl_train),file=paste0("../../data/ML_out/",classif,"/samples_train.txt"), ncolumns = 1)
save(copd_ctrl_train, file = paste("../../data/ML_out/",classif,"/copd_ctrl_train.Rda", sep = ""))

# Train set for mrmr
copd_ctrl_train_mrmr <- copd_ctrl_train
copd_ctrl_train_mrmr$dis_condition <- as.integer(copd_ctrl_train$dis_condition)
write.csv(copd_ctrl_train_mrmr, paste0("../../data/ML_out/",classif,"/mRMR/copd_ctrl_train.csv"), row.names = FALSE)

# Train
print("Class imbalane: train set")
round(prop.table(table(copd_ctrl_train$dis_condition))*100,3)

# Test
copd_ctrl_test <- testing(copd_ctrl_split)
save(copd_ctrl_test, file = paste("../../data/ML_out/",classif,"/copd_ctrl_test.Rda", sep = ""))
print("Class imbalane: test set")
round(prop.table(table(copd_ctrl_test$dis_condition))*100,3)

load(paste("../../data/ML_out/",classif,"/copd_ctrl_test.Rda", sep = ""))
load(paste("../../data/ML_out/",classif,"/copd_ctrl_train.Rda", sep = ""))
#
# write.table(all, "../../data/ML_out/two_class/all.csv",quote=F,sep="\t",row.names=F,col.names=T)
