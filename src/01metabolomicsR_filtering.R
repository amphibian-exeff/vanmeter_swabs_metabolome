# https://xikunhan.github.io/metabolomicsR/docs/articles/Introduction.html

library(metabolomicsR)
library(data.table)
library(ggplot2)
library(cowplot)
library(plotROC)
library(ggstatsplot)
library(M3C) #??

# the metabolite class is presumably already annotated
# Load the metabolomic dataset from an excel file

file_path <- base::system.file("extdata", "QMDiab_metabolomics_OrigScale.xlsx", package = "metabolomicsR", mustWork = TRUE)

df_plasma <- load_excel(path = file_path,
                        data_sheet = 1,
                        feature_sheet = 4,
                        sample_sheet = 8,
                        sampleID = "QMDiab-ID",
                        featureID = "BIOCHEMICAL"
)

df_plasma #S4

# load urine metabolomic data
df_urine <- load_excel(path = file_path,
                       data_sheet = 2,
                       feature_sheet = 5,
                       sample_sheet = 9,
                       sampleID = "QMDiab-ID",
                       featureID = "BIOCHEMICAL"
)

df_urine <- update_Metabolite(df_urine, dataset = "COMP_IDstr", action = "change_featureID")

# load saliva metatabolomic data
df_saliva <- load_excel(path = file_path,
                        data_sheet = 3,
                        feature_sheet = 6,
                        sample_sheet = 10,
                        sampleID = "QMDiab-ID",
                        featureID = "BIOCHEMICAL"
)

#sampleData
View(df_plasma@sampleData)

#assayData
View(df_plasma@assayData) # an object of metabolite
dim(df_plasma@assayData)
colnames(df_plasma@assayData)

#featureData, annotated metabolites and unassigned features w ancillary information
View(df_plasma@featureData)
dim(df_plasma@featureData)
colnames(df_plasma@featureData)

# qc
# winsorize detection of outliers
p <- plot_QC(df_plasma)
p$p
df_plasma_QC <- QC_pipeline(df_plasma, replace_outlier_method = "winsorize", impute_method = NULL)

p <- plot_QC(df_urine)
p$p
df_urine_QC <- QC_pipeline(df_urine, replace_outlier_method = "winsorize", impute_method = NULL)

p <- plot_QC(df_saliva)
p$p
df_saliva_QC <- QC_pipeline(df_saliva, replace_outlier_method = "winsorize", impute_method = NULL)

# transformation
# transformation methods: log (natural logarithm), pareto scale, scale, and rank-based inverse normal transformatio
df_plasma_QC <- impute(df_plasma_QC, method = "half-min")
View(df_plasma_QC)
df_plasma_scale <-  transformation(df_plasma_QC, method = "log")
df_plasma_scale <-  transformation(df_plasma_scale, method = "scale")

# if no features were selected, randomly show 16 metabolites
plot_Metabolite(df_plasma_QC, plot = "boxplot", x = "T2D", color ="ETHNICITY", shape = "T2D")

# select three metabolites
df_plasma_QC@featureData
plot_Metabolite(df_plasma_QC, x = "T2D", plot = "boxplot", feature_name = c("cholesterol",  "choline", "creatine"))

# comparisons between groups using `ggbetweenstats`
plot_Metabolite(df_plasma_QC, x = "T2D", plot = "betweenstats",  feature_name = c("cholesterol",  "choline"))

#boxplot after transformation
plot_Metabolite(df_plasma_scale, plot = "boxplot", x = "T2D", color ="ETHNICITY", shape = "T2D")

#Visualization: histogram
plot_Metabolite(df_plasma_scale, plot = "histogram", color = "T2D")

plot_Metabolite(df_plasma_scale, plot = "histogram", color = "ETHNICITY")

#Normalization
#Normalization of metabolomics data is an important step to remove systematic variation, and to 
#facilitate the identification of real biological metabolites. We provide popular normalization methods: 
#batch-based normalization, QC sample-based/nearest QC sample-based normalization, and LOESS normalization. 
#The latter three methods are useful when the measurements of internal quality control samples are provided. 
#Briefly, in the “batch_norm” function, we implemented a batch-based normalization method that the raw values 
#were divided by the median values of samples in each instrument batch to force the median value to one in 
#each batch. “QCmatrix_norm” function will use reference quality control samples to normalize 
#raw metabolite measurement values.

# Take the urine data for example, the batch of the samples are missing. To illustrate the usage of the normalization method, we assume that 120 samples in each batch. 
unique(featureData(df_urine_QC)$PLATFORM)
## [1] "GC/MS"     "LC/MS Neg" "LC/MS Pos"
# add new columns by reference (column name uppercase)
sampleData(df_urine_QC)[, `GC/MS` := rep(1:3, each = 120)[1:359]]
sampleData(df_urine_QC)[, `LC/MS NEG` := rep(1:3, each = 120)[1:359]]
sampleData(df_urine_QC)[, `LC/MS POS` := rep(1:3, each = 120)[1:359]]

# we also added an injection order column (in the current order)
df_urine_QC <- update_Metabolite(df_urine_QC, sampleData(df_urine_QC)[, 1], action = "injection_order")
## Creat a new column `ID_injection_order`.
# and set one QC sample for each 10 samples. 

sampleData(df_urine_QC)[, QCsample := ifelse(1:359%%10 == 0, 1, 0)]
sampleData(df_urine_QC)[, QCsample := factor(QCsample, levels = c(1, 0))]

sampleData(df_urine_QC)[, sampleID := ifelse(1:359 %% 10 == 0, paste0("MTRX-", sampleID), sampleID)]

table(sampleData(df_urine_QC)$QCsample)
v_features <- featureData(df_urine_QC)$featureID[1:16]
plot_injection_order(df_urine_QC, color = "QCsample", shape = "GC/MS", feature_name = v_features)

assayData(df_urine_QC)$sampleID <- sampleData(df_urine_QC)$sampleID # also change the ID in assayData
df_urine_QC

#We then selected the first 16 metabolites to show boxplot

v_features <- featureData(df_urine_QC)$featureID[1:16]

plot_injection_order(df_urine_QC, color = "QCsample", shape = "GC/MS", feature_name = v_features)

#Normalization 1: batch-norm
df_urine_QC_norm <- batch_norm(df_urine_QC)
## 
## Platform information in @featureData:
## platforms
##     GC/MS LC/MS Neg LC/MS Pos 
##       153       355       278 
## 
## 
## Sample size in platform GC/MSGC/MS
##   1   2   3 
## 120 120 119
p1 <- plot_injection_order(df_urine_QC_norm, color = "QCsample", shape = "GC/MS", feature_name = v_features)
p1

# other normalization methods available...

#Dimensional reduction
#Dimensional reduction strategies on metabolites data can be used to detect batch effects, 
#sample outliers, and real biological subgroups. We included principal components analysis (PCA), manifold 
#approximation and projection (UMAP), and t-distributed stochastic neighbor embedding (tSNE) methods.
df_plasma_PCA <- run_PCA(df_plasma_QC)

plot_PCA(df_plasma_PCA, color ="ETHNICITY", shape = "T2D")
plot_tsne(df_plasma_QC, color ="ETHNICITY", shape = "T2D")

#Correlation
# pairwise correlation of metabolites between two different fluids

# plasma vs urine
dd <- correlation(df_plasma_QC, df_urine_QC, method = "spearman")
## Identify  356  samples in data A, and  359  samples in data B, with  312  overlap samples.
## Identify  545  features in data A, and  786  features in data B, with  248  overlap features.
dd <- merge(dd, featureData(df_plasma_QC), by.x = "term", by.y = "featureID")

p <- ggplot(dd, aes(x = SUPER_PATHWAY, y = r, fill = SUPER_PATHWAY)) +
  geom_boxplot() +
  geom_jitter( size=0.4, alpha=0.9) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set1")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  labs(title = "plasma vs urine")

p
