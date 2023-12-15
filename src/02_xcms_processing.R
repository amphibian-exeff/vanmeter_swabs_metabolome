# Peak Detection 
# centWave algorithm perform peak density and wavelet based chromatographic 
# peak detection for high resolution LC/MS data in centroid m

cwp <- CentWaveParam(ppm = 25, peakwidth = c(5, 10), snthresh = 10) 
xcms_raw_data <- findChromPeaks(swab_metabolites, param = cwp)

# Retention Time Correction You may need to adjust retention times across samples for more reliable comparison.
adjusted_data <- adjustRtime(xcms_raw_data, param = ObiwarpParam(binSize = 0.1))

# Feature Definition Group features that correspond to the same metabolite across samples.
grouped_features <- groupChromPeaks(adjusted_data)

# Load Your Database Load the database that contains known m/z ratios and retention times.
database <- read.csv("path_to_your_database.csv")

# Annotate your detected peaks using the database.
annotated_data <- annotate(grouped_features, database, mzabs = 0.01, mzppm = 15, rt = 5, ppm = 20)

# Review Results Check the annotated peaks.
head(annotated_data@featureData)

#  You will need to replace "path_to_your_data_files" with the actual path to your mass spectrometry data files and "path_to_your_database.csv" with the path to your CSV database.
# The parameters for peak detection and annotation (CentWaveParam, ObiwarpParam, and annotate) need to be chosen based on the specifics of your data and database.
This code is quite generalized and assumes a certain data format. You will need to adapt it to the specifics of your experimental setup and database schema.
The annotate function as used here is illustrative. The actual function you use will depend on the format of your database and how you wish to match the experimental data to the database entries.

Each step's parameters and functions need careful adjustment to match your data's characteristics, and the above is just a very broad outline. Always refer to the xcms documentation and vignettes for detailed examples and explanations of the various functions and their parameters.

Where can I get annotation databases for xcms that contain m/z ratios and retention times?

Annotation databases containing m/z ratios and retention times for use with xcms are crucial for metabolomics studies. There are several publicly available databases that you can use:
  
  METLIN: A large metabolite database that includes m/z values, small molecule and peptide data. It has some retention time information, mainly for LC-MS/MS data.

MassBank: A high-quality mass spectral database that offers spectral data from various instruments and includes some retention indices.

HMDB (Human Metabolome Database): Contains detailed information about small molecule metabolites found in the human body, including chemical, clinical, and molecular biology/biochemistry data along with mass spectrometry data.

LipidMaps: A comprehensive resource for Lipid structures and their MS/MS data, though retention times are not the primary focus.

GNPS (Global Natural Products Social Molecular Networking): A platform for sharing and annotating tandem mass spectrometry data. While it's more focused on molecular networking, it has extensive libraries that may include retention time data.

PubChem: Provides mass spectral information and can be a source for compound identification by m/z ratios. Retention times are not the main data available, but sometimes they are included.

MoNA (MassBank of North America): Similar to MassBank, providing a large repository of mass spectra including retention times when available.

MZmine: While not a database itself, MZmine is an open-source software for mass-spectrometry data processing which can be used alongside databases for metabolite identification.

NIST (National Institute of Standards and Technology): Their mass spectral library includes retention indices and is widely used, although access to the full library may require purchase.

To use these databases with xcms for annotation, you will often need to download the database in a format compatible with R (like CSV or SQLite), or access the database through an API if one is available.
Some databases require you to register and agree to their terms of use before you can download data. After obtaining the data, you might need to process it into a format that matches the requirements of your xcms workflow for peak annotation.

After I have annotated feature data how can I visualize that data and find significant peaks ChatGPT

After annotating features in your mass spectrometry data using xcms, you'll want to visualize and interpret your results, particularly to identify significant peaks that may correspond to compounds of interest. Here's a general approach to visualizing and identifying significant peaks using R:

Plotting the Raw Data: You can start by plotting the raw data to see the overall quality and distribution of peaks.

R

plotChromatogram(raw_data, type = "XIC")

Plotting the Detected Peaks: After peak detection, you can plot the detected peaks to visualize how they overlay with your raw data.

R

plot(raw_data, type = "p", showPeaks = TRUE)

Feature Plot: You can use feature plots to visualize the features that have been grouped across your samples.

R

plotChromPeaks(grouped_features)

EIC (Extracted Ion Chromatogram) Plot: To visualize specific ion traces, you can plot the EIC of annotated features.

R

plotEIC(raw_data, mz = "m/z_value", rtrange = c("start_rt", "end_rt"))

Heatmaps: Heatmaps can be useful to display the intensity of features across different samples.

R

plotHeatmap(grouped_features, method = "mean")

Statistical Analysis for Significant Peaks: To find significant peaks, you can perform statistical analysis. A common approach is to use univariate statistics like t-tests or ANOVA for comparing groups, if you have labeled groups within your data.

R

# Assuming 'grouped_features' contains peak intensity data and 'sample_group' is a factor indicating group membership res <- t.test(grouped_features$intensity ~ sample_group)

Multivariate Analysis: You can also apply multivariate statistics, such as principal component analysis (PCA) or partial least squares discriminant analysis (PLS-DA), to identify peaks that contribute to the separation between groups.

R

library(ropls)
pca_res <- opls(grouped_features$data_matrix, predI = sample_group, plotL = TRUE)

Volcano Plots: For datasets with groups, volcano plots are a great way to visualize both the fold change and statistical significance of features.

R

# Assuming 'stat_results' contains log2 fold changes and p-values with(stat_results, plot(log2FoldChange, -log10(p.value), pch = 20, main = "Volcano Plot"))

Pathway Analysis: If you have identified significant metabolites, you can further analyze their biological relevance by mapping them onto metabolic pathways.

R

# This would typically be done through an external database or web-based tool like MetaboAnalyst.

Please keep in mind that these code snippets are quite generalized and will need to be adapted to fit the structure of your data and the specifics of your experimental design. It's also important to perform appropriate data normalization and transformation before applying statistical tests, to ensure that the results are reliable and meaningful. The packages ggplot2, pheatmap, or other specific packages may offer more sophisticated and customizable plotting options.

Is xcms only for gc mass specs?
  ChatGPT

No, xcms is not only for GC-MS (gas chromatography-mass spectrometry) data; it is a versatile package that can handle data from a variety of mass spectrometry platforms, including LC-MS (liquid chromatography-mass spectrometry), CE-MS (capillary electrophoresis-mass spectrometry), and direct-infusion mass spectrometry.

xcms was originally developed with a focus on metabolomics data, which can be generated from any of these platforms. It provides tools for preprocessing, such as peak detection, retention time alignment, and isotope identification, which are applicable to data from any type of mass spectrometry analysis. This makes it a popular choice for a wide range of applications in both untargeted and targeted metabolomics.

If I have a database of samples with m/z ratios and retention time, what code would I use to annotate metabolites using CAMERA?
  ChatGPT

CAMERA is an R package that complements xcms for the annotation of LC/MS metabolomics data. It helps in the annotation of isotopes, adducts, and fragments in complex mass spectrometry datasets. After processing your data with xcms, you can use CAMERA to further annotate your feature data.

Here is an outline of steps and corresponding R code to annotate metabolites using CAMERA after you've processed your data with xcms:

Load the necessary libraries:

R

library(xcms)
library(CAMERA)

Process your data with xcms: This involves reading the data, peak detection, and retention time correction. Assume you have already done this and have an xcmsSet object named xs.

R

# xs <- your processed xcmsSet object

Annotate with CAMERA: Use the xs object to find isotopes and adducts.

R

an <- xsAnnotate(xs)
an <- groupFWHM(an)
an <- findIsotopes(an)
an <- findAdducts(an)

Compare CAMERA annotations with your database: Suppose your database of m/z ratios and retention times is loaded into a data frame named myDB.

R

myDB <- read.csv("path_to_your_database.csv")

# Assuming your database has columns 'mz' for m/z ratio and 'rt' for retention time for (i in seq_along(an@isotopes)) {
  peak <- an@isotopes[[i]]
  matched <- myDB[abs(myDB$mz - peak$mz) < mzTolerance & abs(myDB$rt -
peak$rt) < rtTolerance, ]
  if (nrow(matched) > 0) {
    # Do something with the matched annotations, like appending them to the peak data
  }
}

In the above pseudo-code, mzTolerance and rtTolerance are the thresholds for matching the m/z and retention times between your annotated peaks and the database. You will need to define these according to the resolution and accuracy of your MS data.

Review and export the results: Once you have matched your CAMERA annotations with your database, you can review the annotations and export the results to a CSV or other file format.

R

# To view the annotation results
summary(an)

# Exporting annotated peak list
write.csv(peakTable(an), "annotated_peak_list.csv")

Please note that this code is a general outline and may need to be adapted for your specific data and the format of your database. The CAMERA package is quite powerful but also complex, so detailed reading of the package documentation and vignettes is recommended to effectively utilize its full functionality.
