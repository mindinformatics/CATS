# make GSEA txt files for all the expression 
source("~/Desktop/CATS/GEO2R/GSEA/expresionfile.R")

library(Biobase)
library(GEOquery)

setwd("~/Desktop/CATS/GEO2R/DatasetCollection/GSE25926_Muller/")
if(!dir.exists("GSEA")){dir.create("GSEA")}

###########################################################
#Excerpt of Astrocytes.R 
#Cases 10: astrocytes p1, p17, p17g, p30
#Controls 28: Astroglia, Neurons, OLS,OPCS, MOGS
gset <- getGEO("GSE9566", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples

gsms <- "00000000010101000000000101100011000011"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

##################################################################