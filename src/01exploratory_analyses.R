#View(swab_metabolites)

# need to bring along the design variables for the analysis

# Preserve existing row names as metabolite/feature IDs
X <- as.matrix(swab_metabolites)

# Confirm its structure
dim(X)          # expected: 8721 rows x 274 columns
head(rownames(X)) # intensities
head(colnames(X)) # sampels

# 274 samples x 8721 metabolite features
X_samples <- t(X)

# After transposition:
rownames(X_samples) #  sample IDs
colnames(X_samples) # metabolite IDs

#Mass-spec intensities are often highly right-skewed, so a log transform is generally useful for exploratory visualizations. 
# log1p() computes log⁡(1+x)log(1+x), which avoids producing -Inf for zero intensities.
X_log <- log1p(X_samples)

# Drop features with no variation across samples (no standard deviation)
keep <- apply(X_log, 2, sd, na.rm = TRUE) > 0
X_log <- X_log[, keep, drop = FALSE]

sample_distance <- dist(X_log)
sample_cluster <- hclust(sample_distance, method = "complete")

plot(sample_cluster, main = "Hierarchical clustering of samples", xlab = "")


# If your matrix has missing values, PCA and many heatmap functions will need them handled first. 
# A simple exploratory approach is to replace missing values in each feature with its observed minimum:
sum(is.na(X_samples))
# there are no NAs

# PCA of samples
pca <- prcomp(X_log, center = TRUE, scale. = TRUE)

variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

plot(
  pca$x[, 1],
  pca$x[, 2],
  pch = 19,
  xlab = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
  ylab = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
  main = "PCA of swab-metabolite samples"
)

head(rownames(pca$x))

# Heatmap of variable features
# Rather than heatmapping all 8,721 features, select the most variable ones across samples:
feature_variance <- apply(X_log, 2, var, na.rm = TRUE)
top_features <- names(sort(feature_variance, decreasing = TRUE))[1:50]
heatmap_matrix <- t(X_log[, top_features, drop = FALSE])

library(pheatmap)

pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = FALSE,
  main = "50 most variable metabolite features"
)

# One boxplot per sample
boxplot(
  log1p(X),
  outline = FALSE,
  las = 2,
  ylab = "log1p intensity",
  main = "Intensity distributions by sample"
)

library(ComplexHeatmap)

Heatmap(
  heatmap_matrix,
  name = "log1p intensity",
  row_title = "Top variable features",
  column_title = "Samples",
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)


# Total signal across metabolite features, per sample
total_signal <- colSums(X, na.rm = TRUE)

plot(
  total_signal,
  type = "b",
  pch = 19,
  xlab = "Sample number",
  ylab = "Total intensity",
  main = "Total metabolite signal by sample"
)


# Fraction of missing values for every sample
missing_fraction <- colMeans(is.na(X))

barplot(
  missing_fraction,
  las = 2,
  ylab = "Fraction missing",
  main = "Missing metabolite intensities by sample"
)
