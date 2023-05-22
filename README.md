# Machine Learning for the characterization of COPD

Master Thesis GitHub repostory (Machine Learning for the characterization of COPD).

## Table of contents
- [Abstract](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main#abstract)
- [Folder Structure](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main#folder-structure)
   - [data](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main/data)
   - [src](https://github.com/airipl/Pose-Lagoa_et_al_COPD_prediction/tree/main/src)
- [Code]()

## Abstract
Chronic Obstructive Pulmonary Disease (COPD) is a complex, heterogeneous, highly prevalent, 
and yet underdiagnosed disease with poor outcomes due to the difficulties of an early diagnosis. 
This study aims to enhance the binary patient classification of COPD using gene expression data 
from the Lung Tissue Research Consortium. To achieve this, we employ various criteria, including 
intrinsic data characteristics (data-driven) and an external information source (DisGeNET),
to identify the most relevant genes. Subsequently, we evaluate the performance of different classifiers, 
such as Random Forest (RF), Support Vector Machines - polynomial and radial kernel -(SVM-poly, SVM-rad), 
k-Nearest Neighbors (kNN), Generalized Linear Models (GLM), and XGBoost (XGB). Our results show that the
data-driven gene selection approach yields the highest accuracies. Specifically, SVM-poly (85,8%) and k
NN (81,6%) are the models that achieve better accuracies in cross-validation and independent test data, 
respectively. Notably, the classification performance is less optimal for controls and COPD patients with 
mild to moderate disease stages.

## Folder Structure

The project has the following folder structure:

``` 
├── data/
│   │  
|   ├── COPDDisgenet.txt    <- COPD-related genes extracted from DisGeNET
|   ├── data_driven.txt    <- data-driven (DEA union mMRM) list of genes
|   ├── eset_genes.txt    <- list of genes of my expression data
|   ├── expansion_data_driven.txt   <- expansion of data-driven (list of genes) 
|   ├── expansion_dgn_genes.txt   <- expansion of copd-related (list of genes) 
|   ├── expansion_intersection.txt   <- intersection of both expansions (list of genes) 
|   ├── expansion_union.txt   <- union of both expansions (list of genes)   
|   ├── minmax_100.txt   <- top 100 genes from mRMR algorithm
|   └──	samples_train.txt   <- sample identifiers of the training set
│
└── src/   
    │
    ├── data/  
    |   |
    |   ├── load_files_GSE47460.R   <- load raw data
    |   └──	omnipath.R    <- extraction of the different collection of genes and inspection of the relations among them
    │
    ├── expression_objects/   
    │   │  
    |   └──	eset_objects_GSE47460.R   <- preprocessing steps for the generation of the Expression Set Objects 
    │
    ├── models_out/         
    │   │                
    |   ├── dea.R   <- Differential Expression Analysis
    |   ├── gsea.R    <- Gene Set Enrichment Analysis
    |   ├── run_ML_mRMR_gs_cv_rep.R   <- ML algorithms using grid search tuning methodoly with a simple 10-fold cross-validation
    |   ├── run_ML_mRMR_gs_cv_simple.R    <- ML algorithms using grid search tuning methodoly with a repeated 10-fold cross-validation
    |   ├── run_ML_mRMR_optm.R    <- ML algorithms using Bayes optimization tuning methodoly with a repeated 10-fold cross-validation
    |   └──	train_test.R    <- generation of train test sets
    │
    ├── random/
    │   │  
    |   ├── outliers.R    <- obtain list of outlier samples
    |   ├── outliers.txt    <- file to run outliers.R
    │   └── random.R    <- randomly generate ML models with using a determined number of genes as input
    │
    ├── support/
    │   │  
    │   └── filter.py    <- filter mRMR result files 	
    │
    ├── eda/
    │   │  
    │   ├── eda_both_plt.R    <- explorative data analysis of expression data of both platforms
    │   └── initial_eda.R    <- explorative data analysis of expression data by platform
    │
    │
    └── visualization_eda/ 
        │  
        ├── RadarChart.R    <- Radar Chart of ML resutls generation
        ├── conclusion.R    <- results comparison between the different ML algorithms
        ├── eda_both_plt.R    <- explorative data analysis of expression data of both platforms
        ├── initial_eda.R    <- explorative data analysis of expression data by platform
        ├── metrics.R    <- recolect results (metrics, missclasified samples) of ML models
        ├── metrics_comparison.R    <- comparison between the different ML algorithms
        ├── miss_sample.R    <- analysis of missclasified samples
        ├── pubmed_publications.R    <- obtain PubMed articles associated with gene-COPD terms
        ├── random_analysis.R    <- analysis of random ML models results
        └── time_analysis.R     <- time analysis of different ML models
```

## Code
### Data loading and preparation
First, we downloaded the raw data from GEO and we prepared it for reading and processing (`src/data/load_files_GSE47460.R`). We also applied a pipeline for data preprocessing (normalization, bg correction, filtering of genes and probes...) and generation of the Expression Set objects (`src/expression_objects/eset_objects_GSE47460.R`). The script `outliers.R` was generated for the detection of the outlier samples. A first explorative data analysis was computed to each platform separated (`src/eda/initial_eda.R`) and joining both of them (`src/eda/eda_both_plt.R`) (both scripts update the expression objects if necessary).
