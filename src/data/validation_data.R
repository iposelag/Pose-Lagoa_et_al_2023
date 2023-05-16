###################################
## Load data for validation task ##
###################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required packages
library(sva)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data

# Metadata
meta_data <- read.table("../../data/raw/IriaPose/MetaData.txt", header = TRUE)
# Gene expression
gene_expression <- read.table("../../data/raw/IriaPose/GeneExpression.txt", header = TRUE)
# Check that the samples are the same in both data sets
table(colnames(gene_expression) %in% meta_data$Sample)
# Switch - by . and check that the samples are the same
meta_data$Sample <- gsub("-", "\\.", meta_data$Sample)
table(colnames(gene_expression) == meta_data$Sample)

## ----------------------------------------------------------------------------------------------------------------------------------------------
# BATCH effect correction (platform): ComBAT
batch = meta_data$Study
modcombat = model.matrix(~ Category, data=meta_data)
combat_edata = ComBat(dat=gene_expression, batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)

# Data frame of interest with ComBat data
combat_edata_ph <- as.data.frame(t(combat_edata))
dim(combat_edata_ph)

combat_edata_ph$dis_condition <- as.factor(meta_data$Category)
combat_copd_ctrl <- combat_edata_ph
dim(combat_copd_ctrl)
combat_copd_ctrl = combat_copd_ctrl[ ,
                                     c("dis_condition",names(combat_copd_ctrl)[
                                       names(combat_copd_ctrl) != "dis_condition"])]
dim(combat_copd_ctrl)

# Change levels factor names
combat_copd_ctrl$dis_condition
new_levels <- c("CTRL", "COPD")
# Check that are correct changed
table(combat_copd_ctrl$dis_condition)

# Cambiar los niveles de la variable factor
levels(combat_copd_ctrl$dis_condition) <- new_levels
# Check that are correct changed
table(combat_copd_ctrl$dis_condition)

combat_copd_ctrl$dis_condition <- relevel(combat_copd_ctrl$dis_condition,
                                          ref="COPD")
# write.table(combat_copd_ctrl, "../../data/ML_out/validation_data/copd_ctrl_jon.csv",quote=F,sep="\t",row.names=F,col.names=T)
save(combat_copd_ctrl, file = "../../data/ML_out/validation_data/copd_ctrl_jon.Rda")
