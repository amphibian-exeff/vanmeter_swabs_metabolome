# Van Meter salamander swab metabolome data



#Install and load supporting libraries.
print(Sys.info()[4])
#BiocManager::install("pathview")

library(dplyr)
library(forcats)
library(FSA)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(BGGM)
library(ggm)
library(corrplot)
library(qgraph)
library(robFitConGraph)
library(reshape2)
library(glasso)
library(igraph)
library(GGally)
library(matrixcalc)
library(Matrix)
library(pathview)
library(pheatmap)
remotes::install_github("XikunHan/metabolomicsR")
library(metabolomicsR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("xcms")

#browseVignettes("pathview")

print("list of loaded packages: ")
print((.packages()))

#tom epa windows
if(Sys.info()[4]=="DZ2626UTPURUCKE"){
  rvm_root <- file.path("c:", "git", "vanmeter_swabs_metabolome")
}
if(Sys.info()[4]=="LZ2626UTPURUCKE"){
  rvm_root <- file.path("c:","git","vanmeter_swabs_metabolome")
}
if(Sys.info()[4]=="LZ26TPURUCKE-2"){ 
  # tom windows 2023 laptop
  rvm_root <- file.path("c:", "Users", "tpurucke", "git", "vanmeter_swabs_metabolome")
}

print(paste("Root directory location: ", rvm_root, sep=""))

rvm_data_in <- file.path(rvm_root, "data_in")
rvm_data_out <- file.path(rvm_root, "data_out")
rvm_graphics <- file.path(rvm_root, "graphics")

#check to see if directories are accessible
boo = file.exists(file.path(rvm_data_in,"/vanmeter_swabs_metabolome.csv"))
print(paste("check to see if R can access GSF file OK: ", boo))


rvm_swabs_t <- read.csv(file.path(rvm_data_in,"/vanmeter_swabs_metabolome.csv"), stringsAsFactors = TRUE)
dim(rvm_swabs_t)

rvm_swabs_temp <- t(rvm_swabs_t)
dim(rvm_swabs_temp)
#View(rvm_swabs_temp)
colnames(rvm_swabs_temp) <- rvm_swabs_temp[1,]


rvm_swabs <- rvm_swabs_temp[-1,]
rvm_swabs <- cbind (as.vector(rownames(rvm_swabs)),rvm_swabs)
colnames(rvm_swabs)[1] <- "sample_id"
colnames(rvm_swabs)[1:7]
dim(rvm_swabs)
