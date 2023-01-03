# GTEx_SmokerStatus - Brain Data

###About
Smoking a commonly known confounder for many diseases can affect gene expression and regulation in several prtions of the brain. Being closely linked to the nicotine reward pathway, excessive smoking is known to increase activity of genes that code for proteins responsible for activity of Neuronal Nicotinic Acetycholine Receptors (nAChRs). This change in gene expression is usually irreversible and allows us to understand differences in brain function over time due to smoking. These signatures remain in individuals and their smoking status can be predicted from this. Exploratory Analysis of changes in gene expression in the Brain based on transcripts per million and covariates data obtained from GTEx Portal could help us test existing literature and add information about confounding principles like smoking. The portal provides gene expression information for different sub sections of the brain. The following brain regions data is analysed :

1. Substantia Nigra
2. Hypothalamus
3. Hippocampus
4. Nucleus Accumbens
5. Pre-frontal cortex
6. Amygdala


###Data sets used :
Every above brain sub-region has a unique transcripts per million gene expression file (tpm file) and a covariates matrix associated with it. A gene regualtion matrix  ("Regulation of Gene expression.csv") has information about the expected changes in gene expression for the choosen set of genes in each brain sub-region on different levels of exposure to smoking. The tpm file that shows gene expression may not capture smoking alone as a signature and could be affected by the several other known and unknown variables. The covariates however could be a better indicator for capturing differences in gene expression with respect to smoking alone. 

###Code files:
Majority of the R code files have the naming pattern of "brain_region_clustering" for example "Amygdala_clustering". These files contain the R scripts used to cluster data using K-means algorithm and generate plots of the clusters made. Moreover, it has the code used to understand the differences between the clusters generated. 

The "corr.R" file is a script to generate correlation matrices for each brain sub-region. The correlation matrix has the degree of correlation between each covariate from the covariates matrix and the genes expression from the tpm file of the genes of interest choosen from the gene regulation matrix for each brain region.

The pdf files have plots of the clusters generated. They name indicates the brain region and the number of covariates used to make the clusters. The covariates with high correlation with the genes of interest are choosen for this clustering. 

The output csv files associated with each brain region have results of cluster comparisons for each clustering model (when different number of covariates are used, slightly different clustering results are obtained). 
