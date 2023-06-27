
input <- c("data_driven", "data_driven_expansion", "disgenet", "disgenet_expansion", 
           "expansion_intersection", "expansion_union")
for(i in input){
  # Load miss_samples table
  miss_samples <- read.csv(paste("../../data/ML_out/two_class/optm/",i,
                                    "/miss_samples/miss_samples_train.csv", sep=""))
  miss_samples <- setNames(miss_samples, c(paste0("rf_",i),
                                           "sample_id",
                                           paste0("svm_r_",i),
                                           paste0("svm_p_",i),
                                           paste0("knn_",i),
                                           paste0("glm_",i),
                                           paste0("xgb_",i)))
  assign(paste0("miss_samples_",i), miss_samples)
}

# Crea una lista con todos los data frames
df_list <- list(miss_samples_data_driven, miss_samples_data_driven_expansion, 
                miss_samples_disgenet, miss_samples_disgenet_expansion, 
                miss_samples_expansion_intersection, miss_samples_expansion_union)

# Une los data frames por la columna "sample_id"
merged_df <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE, sort = FALSE), df_list)

# Save the data frame
write.csv(merged_df, paste("../../data/ML_out/two_class/optm/tables/miss_samples_results.csv",sep=""), 
          row.names = FALSE)

## Load phenotypic characteristics
load("../../data/ExpressionSetObjects/eset.Rda")
dim(eset)
n_samples <- ncol(eset)
n_genes <- nrow(eset)
phenotypic_ctics <- pData(eset)

interest <- c("dis_condition","GOLD_stage", "sex", "smoker", "age", 
              "pneumocystis_colonization", "emphysema", "pred_dlco", 
              "pred_fev_prebd", "pred_fev_postbd", "pred_fvc_prebd", 
              "pred_fvc_postbd", "platform_id")
phenotypic_ctics_interest <- phenotypic_ctics[, interest]
phenotypic_ctics_interest <- phenotypic_ctics_interest %>%
  mutate(sex = factor(sex),
         pneumocystis_colonization = factor(pneumocystis_colonization),
         smoker = factor(smoker),
         GOLD_stage = factor(GOLD_stage),
         platform_id = factor(platform_id),
         age = as.numeric(age),
         pred_dlco = as.numeric(pred_dlco),
         pred_fev_prebd = as.numeric(pred_fev_prebd),
         pred_fev_postbd = as.numeric(pred_fev_postbd),
         pred_fvc_prebd = as.numeric(pred_fvc_prebd),
         pred_fvc_postbd = as.numeric(pred_fvc_postbd))
phenotypic_ctics_interest$sample_id <- rownames(phenotypic_ctics_interest)

# Save all train (test) miss classified samples
copd_ctrl_condition <- subset(phenotypic_ctics_interest, sample_id %in%
                                merged_df$sample_id)
merged_df_pheno <- merge(merged_df, copd_ctrl_condition, by="sample_id")

# Save the data frame
write.csv(merged_df_pheno, paste("../../data/ML_out/two_class/optm/tables/miss_samples_results_pheno.csv",sep=""), 
          row.names = FALSE)
