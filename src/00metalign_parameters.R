# metalign
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02159-0

# MetAlign is a software tool developed for the preprocessing and comparison of full scan nominal 
# mass LC-MS data. It was designed to handle large sets of high-resolution MS data, such as the 
# kind obtained from metabolomics experiments.

# In metabolomics, researchers analyze the small molecules (metabolites) present in biological 
# samples. Given the complexity and dynamic nature of the metabolome, preprocessing and analysis 
# of the data are crucial steps. MetAlign's main functions in this context include
# Noise Estimation and Removal: It can distinguish and remove noise from the actual signal, which 
# improves the quality of the resulting data.
# Alignment of Data: Over different runs or samples, the retention time for a particular compound 
# might slightly vary. MetAlign aligns these variations, making it easier to compare samples.
# Data Reduction: It reduces the data size without losing the essential information, which 
# makes subsequent data processing and analysis more manageable.
# Peak Detection: MetAlign identifies and quantifies metabolite peaks in the data.

# Settings used in metalign, per Donna:
# So far the data has gone through metalign only (baseline correction, scaling, and normalization) 
# to get the values in the CSV file

## Accurate or nominal form?
# Data type: accurate mass data
# Mass resolution: 8000
# Amplitude range: 100-10000000
# Echo suppression: Interval 0.45, Percentage 5, 
# Forest suppression: Interal around mass peak 6, Percentage of amplitude of mass peak 3, interval offset 0.5

## Program correction, dataset selection, and baseline correction





# Once preprocessed using MetAlign, the data can be fed into other software tools or statistical 
# programs for further analysis, such as principal component analysis (PCA) or multivariate 
# statistical analysis.

# After running MetAlign for preprocessing of metabolomics data, subsequent data processing steps 
# aim to extract meaningful biological information from the cleaned and aligned data. Here's a 
# general workflow of what's typically implemented:
# Feature Detection and Quantification: If not done during preprocessing, detect peaks or features 
# in the aligned dataset. Quantify the intensity of these features across different samples.
# Normalization: To correct for variations in sample amounts or instrument sensitivity, data is 
# often normalized. Common methods include total sum scaling, normalization to an internal standard, 
# or probabilistic quotient normalization.
# Identification of Metabolites: Match the detected features (based on mass-to-charge ratio and 
# sometimes retention time) to known metabolites in databases such as METLIN, HMDB, or MassBank.Further 
# validation might be required using MS/MS data or standards.
# Statistical Analysis: Univariate statistics, such as t-tests or ANOVA, to identify significantly 
# altered metabolites between conditions.
# Multivariate statistics, like principal component analysis (PCA) or partial least squares-discriminant 
# analysis (PLS-DA), to assess the overall variation in the dataset and group differences.
# Biomarker Discovery: If the study's aim is to find potential biomarkers, features that discriminate 
# between conditions or stages of a disease are further investigated.
# Pathway Analysis: Using identified metabolites, pathway enrichment analyses can be performed to 
# determine which metabolic pathways are affected. Tools like MetaboAnalyst or KEGG can assist in this.
# Data Integration: If data from other omics studies (e.g., transcriptomics, proteomics) is available, 
# integrative analyses can be performed to gain a holistic understanding of the biological system.
# Validation: Important findings, especially potential biomarkers, are often validated using an 
# independent dataset or different analytical techniques.

# Throughout these steps, visualization tools and techniques play a crucial role. Heatmaps, 
# volcano plots, PCA score plots, and box plots are some of the commonly used visual representations 
# in metabolomics data analysis.
