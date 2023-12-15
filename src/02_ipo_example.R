# IPO
# https://www.bioconductor.org/packages/devel/bioc/vignettes/IPO/inst/doc/IPO.html

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("IPO")

# for examples of peak picking parameter optimization:
BiocManager::install("msdata")
# for examples of optimization of retention time correction and grouping
# parameters:
BiocManager::install("faahKO")

library(IPO)

# raw data
# xcms handles the file processing hence all files can be used that can be processed by xcms.
datapath <- system.file("cdf", package = "faahKO")
datafiles <- list.files(datapath, recursive = TRUE, full.names = TRUE)
                        
peakpickingParameters <- getDefaultXcmsSetStartingParams('matchedFilter')
#setting levels for step to 0.2 and 0.3 (hence 0.25 is the center point)
peakpickingParameters$step <- c(0.2, 0.3)
peakpickingParameters$fwhm <- c(40, 50)
#setting only one value for steps therefore this parameter is not optimized
peakpickingParameters$steps <- 2

time.xcmsSet <- system.time({ # measuring time
  resultPeakpicking <- 
    optimizeXcmsSet(files = datafiles[1:2], 
                    params = peakpickingParameters, 
                    nSlaves = 1, 
                    subdir = NULL,
                    plot = TRUE)
})