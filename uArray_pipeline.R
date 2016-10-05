# Last version of this script written on 25/SEP/2016
# by Yann Lambert, MIND Informatics
# ------------------

setwd("/Users/ylambert/Desktop/AMP-AD/Pipelines/")

library(data.table)
library(dplyr)
library(limma)

rm(list = ls())

outputFolder <- "Results/uArray/"
logFolder <- "Results/uArray/Logs/"
processFolder <- "Scripts/"

log <- "Expression microarray pipeline script"
logPush <- function(...){
    log[length(log)+1] <<- paste0(...)
    message(...)
}

outputFileRegistry <- data.table()

# 1) Get dataset information

datasets <- fread("uArray_datasets.csv")
dsNameList <- datasets$Dataset

# 2) Apply pipeline to each dataset
for (ds in dsNameList) {
    dsProperties <- datasets[Dataset == ds]

    outputFileEntry <- list(
        FileName = "",
        DateCreated = "",
        Dataset = ds,
        Study = dsProperties$Study,
        BrainRegion = dsProperties$BrainRegion,
        Stratification = "",
        Contrast = "",
        DataFile = paste0(dsProperties$DataFolder, dsProperties$DataFileName),
        ProcessScript = dsProperties$ProcessScript,
        NProbes = NA,
        NSubjects = NA,
        NRemovedSubjects = NA,
        SubjectDistribution = ""
    )
    
    logPush("Run on ", Sys.time())
    logPush("Dataset: ", dsProperties$Dataset)
    
    # 3) Load and process data files
    source(paste0(processFolder, dsProperties$ProcessScript))
    expSet <- processDataset(dsProperties)
    
    # For each stratification factor:
    sFactors <- strsplit(dsProperties$Stratification, "/")[[1]]
    for (sf in sFactors) {
        
        logPush("Stratification: ", sf)
        
        # 4) Check processed data
        source(paste0(processFolder, "uArray_data_check.R"))
        mData <- dataCheck(expSet, sf)
        
        # 5) Analysis
        
        # 5.1) Create design matrix
        f <- factor(sapply(colnames(mData), function(x) expSet$strat[[sf]]$data[SubjectID == x][[sf]]))
        design <- model.matrix(~f+0)
        colnames(design) <- expSet$strat[[sf]]$levels
        
        subjectDistribution <- paste0(paste(expSet$strat[[sf]]$levels, table(f), sep = ": "), collapse=", ")
        logPush("Subject distribution: ", subjectDistribution)
        
        # 5.2) Fit model et make contrasts
        fit <- lmFit(mData,design)
        fit.c <- contrasts.fit(fit, makeContrasts(contrasts =  expSet$strat[[sf]]$contrasts, levels=design))
        
        # 5.3) Adjust p-values with BH method
        fit.eb <- eBayes(fit.c, trend = T)

        # 5.4) Gene annotation
        ids <- rownames(fit.eb$p.value)
        annotation <- merge(data.table(ProbeID = ids), expSet$annot, by="ProbeID", all.x = T, sort = F)
        fit.eb$genes  <- annotation
        
        # 5.5) Output
        for (contrast in expSet$strat[[sf]]$contrasts) {
            tT <- topTable(fit.eb, sort.by="B", number=nrow(fit.eb$p.value), coef=contrast)
            outputFile <- paste0(outputFolder, ds, "_", sf, "_", contrast, ".csv")
            write.csv(tT, file=outputFile, row.names=F)
            
            logPush("Output file: ", outputFile)
            
            outputFileEntry$Stratification <- sf
            outputFileEntry$SubjectDistribution <- subjectDistribution
            outputFileEntry$Contrast <- contrast
            outputFileEntry$FileName <- outputFile
            outputFileEntry$DateCreated <- Sys.time()
            
            outputFileRegistry <- rbind(outputFileRegistry, outputFileEntry)
        }
    }
    
    write(log, file=paste0(logFolder, ds, ".log"))
}

write.csv(outputFileRegistry, paste0("uArray-pipeline_OutputRegistry.csv"), row.names = F)

