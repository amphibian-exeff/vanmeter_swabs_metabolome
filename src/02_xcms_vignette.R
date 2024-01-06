# https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms-lcms-ms.html

library(xcms)
BiocManager::install("MsExperiment")
library(MsExperiment)
browseVignettes("MsExperiment")
BiocManager::install("msdata")
# library(msdata)

browseVignettes("xcms")

dda_file <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
                        package = "msdata")
dda_data <- readMsExperiment(dda_file)
chr <- chromatogram(dda_data, aggregationFun = "sum", msLevel = 1L)
