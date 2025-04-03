#Prefiltering metabolomics data is crucial for reducing the complexity and removing irrelevant or low-quality 
#information before in-depth analysis. Here are some common approaches for prefiltering metabolomics data:
#Intensity Thresholding: Remove features with peak intensities below a certain threshold. This helps in eliminating low-intensity noise.
#Variability Filtering: Features with a low coefficient of variation
#(CV) across replicates are retained, while those with high variability (suggesting they might be noise) can be filtered out.
#Quality Control (QC) Sample-Based Filtering: Using QC samples (pooled samples analyzed at regular intervals throughout the run) 
#to assess the consistency of feature detection. Features that are not reliably detected in QC samples can be discarded.
#Retention Time (RT) Window: Metabolites elute at different times in chromatographic metabolomics. By defining a retention 
#time window, one can filter out features that elute too early or too late, which are often associated with noise or artifacts.
#Blank Sample Filtering: Comparing the metabolite profiles of actual samples to those of blank samples (containing no biological material).
#Features present in the blank may arise from the instrumentation or sample preparation process and can be removed.
#Mass Defect Filtering: This exploits the subtle differences in exact masses between naturally occurring isotopes. Some noise or 
#artifacts might not adhere to natural isotopic patterns and can thus be identified and filtered out.
#Signal-to-Noise Ratio: Filtering based on the signal-to-noise ratio can help to retain only features that have a strong signal 
#compared to the noise level.
#Known Contaminant Removal: Some contaminants are known to be introduced during sample preparation or analysis. If these are 
#identified in advance, their corresponding features can be removed.
#Database Matching: Only retaining features that can be matched to known metabolites in databases, though this might lead 
#to the exclusion of novel or unidentified compounds.
#Normalization: While not strictly a filtering step, normalizing data can help in reducing systemic variation between samples, 
#making subsequent filtering more effective.

#It's important to note that the choice of prefiltering methods should be informed by the specifics of the 
#metabolomics experiment, the platform used, and the goals of the analysis. Overzealous filtering might lead to the 
#loss of valuable information, while insufficient filtering can cloud the analysis with irrelevant or noisy data.

#Intensity threshold filtering is a common method to remove low-intensity features, which might be noise or irrelevant signals.
#Here's a basic step-by-step guide on how to perform intensity threshold filtering for metabolomic data:
#1. Understand the Data: Before applying any filtering, ensure you're familiar with the structure and characteristics of your data.
#Typically, metabolomic data consists of a matrix where rows represent individual metabolites or features and columns represent 
#individual samples.

# Settings used in metalign, per Donna:
# So far the data has gone through metalign only (baseline correction, scaling, and normalization)

# transpose the matrix
dim(rvm_swabs)
swab_metabolites <- as.data.frame(t(rvm_swabs[,8:8728]))
dim(swab_metabolites)
summary(swab_metabolites)
swab_metabolites[] <- lapply(swab_metabolites, function(x) as.numeric(as.character(x)))
dim(swab_metabolites)
summary(swab_metabolites)
View(swab_metabolites)

# Noise Estimation and Removal: It can distinguish and remove noise from the actual signal, which 
# improves the quality of the resulting data.
# Determine the Threshold: Decide on a threshold value for the minimum intensity. This could be based on:
# A fixed intensity value, often determined empirically.
# A multiple of the median or mean intensity of all features.
# A certain percentile of the intensity distribution.

max_metabolite_intensities <- apply(swab_metabolites, 1, max)
sum_metabolite_intensities <- apply(swab_metabolites, 1, sum)
min(max_metabolite_intensities)
mean(max_metabolite_intensities)
median(max_metabolite_intensities)
max(max_metabolite_intensities)
hist(log(max_metabolite_intensities), breaks=1000)
hist(log(sum_metabolite_intensities), breaks=1000)

# Apply the Threshold: Filter out rows (metabolites/features) where the maximum intensity across all 
# samples is below the chosen threshold.

# metalign only (baseline correction, scaling, and normalization)

# 1% of max as threshold
threshold_01 <- max(max_metabolite_intensities) * 0.01
sum(max_metabolite_intensities > threshold_01)/length(max_metabolite_intensities) # 0.0323
# too aggressive, drops almost 97% of data


