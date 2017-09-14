# Testing Download Code
source('~/Desktop/MousePipeline/Code/Download.R')

setwd('~/Desktop/MousePipeline/Download/')

# clear all the files
system(command = "./remove.sh")

# all current datasets 
Aydin = DownloadFunction("GSE25926")
Aydin = DownloadFirstElement("GSE25926")


Tsai = DownloadFunction("GSE65159")

########################################################################

# Barres has multipule platforms
Bares = DownloadFunction("GSE9566")

GEOtag = "GSE9566"

#are they the same samples?
plat1 = sampleNames(phenoData(GEOobject$`GSE9566-GPL1261_series_matrix.txt.gz`))
plat2 = sampleNames(phenoData(GEOobject$`GSE9566-GPL6096_series_matrix.txt.gz`))

# no overlap in samles
table(plat1 %in% plat2)
table(plat2 %in% plat1)

########################################################################

Pfizer = DownloadFunction("GSE31624")
BMS_GSE56772 = DownloadFunction("GSE56772")
BMS_GSE57528 = DownloadFunction("GSE57528")
BMS_GSE57583 = DownloadFunction("GSE57583")