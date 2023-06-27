# Machine Learning for the characterization of COPD

Master Thesis GitHub repository (Machine Learning for the characterization of COPD).

## Table of contents
- [Abstract](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main#abstract)
- [Folder Structure](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main#folder-structure)
   - [data](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main/data)
   - [src](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main/src)
- [Source Code (src)](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main#source-code-src)

## Abstract

Chronic Obstructive Pulmonary Disease (COPD) is a complex, heterogeneous, highly prevalent, and yet underdiagnosed disease with poor outcomes due to the difficulties of an early diagnosis. This study aims to enhance the binary patient classification of COPD using gene expression data from the Lung Tissue Research Consortium. To achieve this, we employ various feature selection criteria, including intrinsic data characteristics (data-driven), an external information source (curated COPD-related genes), and their respective biological expansions to identify the most relevant genes. Subsequently, we evaluate the performance of different classifiers: Random Forest (RF), Support Vector Machines - polynomial and radial kernel -(SVM-poly, SVM-rad), k-Nearest Neighbors (kNN), Generalized Linear Models (GLM), and XGBoost (XGB). Our results show that the data-driven and curated COPD-related OmniPath expansion gene selection approaches yield the highest performances in cross-validation and independent test data. Specifically, kNN (84,8\%) and SVM-poly (81,6\%) are the models that achieve better accuracies in cross-validation and independent test data, respectively.

## Folder Structure

The project has the following folder structure:

``` 
├── data/
│   │  
|   ├── Supplementary Data <- lists of genes (used as inputs)
|   └── samples_train.txt   <- sample identifiers of the training set
│
└── src/   
    │
    ├── data/  
    |   |
    |   └── load_files_GSE47460.R   <- load raw data
    │
    ├── eda/
    │   │  
    │   ├── eda_both_plt.R    <- explorative data analysis of expression data of both platforms
    │   └── initial_eda.R    <- explorative data analysis of expression data by platform
    │ 
    ├── expression_objects/   
    │   │  
    |   └──	eset_objects_GSE47460.R   <- preprocessing steps for the generation of the Expression Set Objects 
    │   
    ├── feature_selection/
    │   │
    |   ├── mrmr
    │        ├── bash_mrmr.sh
    │        ├── bash_mrmr_random.sh
    │        ├── filter.py
    |        └── randomize_sample_labels_mrmr.R
    |   ├── dea.R   <- Differential Expression Analysis
    |   ├── gsea.R    <- Gene Set Enrichment Analysis
    |   ├── EnrichmentAnalysisAndPathwayMap.R <- Enrichment analysis of curated COPD-related, DEA and mRMR genes
    |   └──	omnipath.R    <- extraction of the different collection of genes and inspection of the relations among them
    | 
    ├── ML_models/         
    │   │                
    |   ├── run_ML_mRMR_gs_cv_rep.R   <- ML algorithms using grid search tuning methodoly with a simple 10-fold cross-validation
    |   ├── run_ML_mRMR_gs_cv_simple.R    <- ML algorithms using grid search tuning methodoly with a repeated 10-fold cross-validation
    |   ├── run_ML_mRMR_optm.R    <- ML algorithms using Bayes optimization tuning methodoly with a repeated 10-fold cross-validation
    |   └──	train_test.R    <- generation of train test sets
    │
    ├── random/
    │   │  
    |   ├── outliers.R    <- obtain list of outlier samples
    |   ├── run_outliers.txt    <- file to run outliers.R
    │   └── random.R    <- randomly generate ML models with using a determined number of genes as input
    │
    └── ml_results_analysis/ 
        │
        ├── RadarChart_ml_input_comparisons.R    <- RadarChart with metrics results
        ├── all_metrics_results.R    <- prepare data for RadarCHart
        ├── analisis_descriptivo.R    <- basic analysis of misclassified samples
        ├── barplot_ml_input_comparison.R   <- Barplots for the comparison among the different ML methods and using different inputs gene sets
        ├── enrichment.R    <- enrichement analysis of misclassified samples
        ├── metrics.R    <- metrics recopilation of ML models
        ├── miss_sample.R    <- analysis of missclasified samples
        ├── miss_samples_heatmap.R    <- heatmap of missclasified sample
        ├── pubmed_publications.R    <- obtain PubMed articles associated with gene-COPD terms
        ├── random_analysis.R    <- analysis of random ML models results
        ├── time_analysis.R    <- Comparison between time of tuning and training by input gene sets and ML model
        └── tuning_methodologies_comparison.R     <- barplot for the comparison between tuning methodologies
```

## Source Code (src)
### Data loading and preparation
First, we downloaded the raw data from GEO and we prepared it for reading and processing (`src/data/load_files_GSE47460.R`). We also applied a pipeline for data preprocessing (normalization, bg correction, filtering of genes and probes...) and generation of the Expression Set objects (`src/expression_objects/eset_objects_GSE47460.R`). The script `outliers.R` was generated for the detection of the outlier samples. A first explorative data analysis was computed to each platform separated (`src/eda/initial_eda.R`) and joining both of them (`src/eda/eda_both_plt.R`) (both scripts update the expression objects if necessary).

### Feature Selection 
Then, we split the data into training and test set (`src/ML_models/train_test.R`). Over the training data, we performed the Differential Expression Analysis (`src/feature_selection/dea.R`), mRMR algorithm (`src/feature_selection/mrmr/`), and obtained the data_driven, COPD-related and expansions lists of genes (`src/feature_selection/omnipath.R`). Moreover, we performed a Gene Set Enrichment Analysis (`src/feature_selection/gsea.R`) and Enrichemnet Analysis (`src/features_selection/omnipath.R`).

### Machine Learning models
Once we have our input genes selected, we run ML models using a grid search using simple 10-fold cross-validation (`src/ML_models/run_ML_mRMR_gs_cv_simple.R`), repeated 10-fold cross-validation (`src/ML_models/run_ML_mRMR_gs_cv_rep.R`) and Bayes optimization with repeated 10-fold cross-validation (`src/ML_models/run_ML_mRMR_optm.R`) for hyperparameter tuning. For evalutaing the significance of our results we run randomizations of the different models mantaining the size of the different input sets (`src/random/random.R`)

### Results analysis
We analysed the results of our ML models performing different comparisons and plotting the results. First we recopilated the ML metrics on .csv files (`src/ml_results_analysis/metrics.R`, `src/ml_results_analysis/all_metrics.R`). Then we generated different plots of metrics results (`src/ml_results_analysis/tuning_methodologies_comparison.R`,`src/ml_results_analysis/barplot_ml_input_comparison.R`, `src/ml_results_analysis/random_analysis.R`, `src/ml_results_analysis/RadarChart_ml_input_comparisons.R`, `src/ml_results_analysis/time_analysis.R`). We also analyzed the misclassified samples (`src/ml_results_analysis/miss_samples.R`,`src/ml_results_analysis/miss_samples_heatmap.R`). The file `src/ml_results_analysis/pubmed_publications.R` was used analyze the number of publications that cited our seed genes along with COPD in PubMed.
