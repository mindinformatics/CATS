# Last version of this script written on 25/SEP/2016
# by Yann Lambert, MIND Informatics
# ------------------

setwd("/Users/ylambert/Desktop/AMP-AD/Pipelines/")

library(data.table)
library(dplyr)
library(limma)
library(edgeR)

rm(list = ls())

outputFolder <- "Results/RNASeq/"
logFolder <- "Results/RNASeq/Logs/"
processFolder <- "Scripts/"

log <- "RNASeq pipeline script"
logPush <- function(...){
    log[length(log)+1] <<- paste0(...)
    message(...)
}

outputFileRegistry <- data.table()

# 1) Get dataset information

datasets <- fread("RNASeq_datasets.csv")
dsNameList <- datasets$Dataset

annot <- fread("Data/ENSEMBL_GRCh38.p7_GenesOnly.csv")

# 2) Apply pipeline to each dataset
for (ds in dsNameList) {
        dsProperties <- datasets[Dataset == ds]
        dtype <- dsProperties$DataType
        
        outputFileEntry <- list(
            FileName = "",
            DateCreated = "",
            Dataset = ds,
            Study = dsProperties$Study,
            DataType = dsProperties$DataType,
            BrainRegion = dsProperties$BrainRegion,
            Stratification = "",
            Contrast = "",
            DataFile = paste0(dsProperties$DataFolder, dsProperties$DataFileName),
            ProcessScript = dsProperties$ProcessScript,
            NGenes = NA,
            NSubjects = NA,
            NRemovedSubjects = NA,
            SubjectDistribution = ""
        )
        
        logPush("Run on ", Sys.time())
        logPush("Dataset: ", ds)
        logPush("Data type: ", dtype)
        
        # 3) Load and process data files
        source(paste0(processFolder, dsProperties$ProcessScript))
        expSet <- processDataset(dsProperties)
        
        # For each stratification factor:
        sFactors <- strsplit(dsProperties$Stratification, "/")[[1]]
        for (sf in sFactors) {
            
            logPush("Stratification: ", sf)

            # 4) Check processed data
            source(paste0(processFolder, "RNASeq_data_check.R"))
            mData <- dataCheck(expSet, sf)

            # 5) Analysis

            # 5.1) Create design matrix
            f <- factor(sapply(colnames(mData), function(x) expSet$strat[[sf]]$data[SubjectID == x][[sf]]))
            design <- model.matrix(~f+0)
            colnames(design) <- expSet$strat[[sf]]$levels
            subjectDistribution <- paste0(paste(expSet$strat[[sf]]$levels, table(f), sep = ": "), collapse=", ")
            logPush("Subject distribution: ", subjectDistribution)

            # 5.2) Transform data according to data type
            if (dtype == "RC") {
                # TMM Normalization with EdgeR package
                dge <- DGEList(counts=mData)
                dge <- calcNormFactors(dge)
                mData <- voom(dge,design,plot=F)
            }
            else if (dtype == "CPM") {
                # Pre-TMM normalized data
                mData <- voom(mData,design,plot=F)
            }
            else if (dtype == "FPKM") {
                # FPKM values must not be used with voom
                # Instead, add a small intercept and use log2
                # https://support.bioconductor.org/p/56275/
                mData <- log2(mData + 0.05)
            }

            # 5.3) Fit model et make contrasts
            fit <- lmFit(mData,design)
            fit.c <- contrasts.fit(fit, makeContrasts(contrasts =  expSet$strat[[sf]]$contrasts, levels=design))

            # 5.4) Adjust p-values with BH method
            fit.eb <- eBayes(fit.c, trend = T)

            # 5.5) Annotation
            geneIDs <- rownames(fit.eb$p.value)
            annotation <- merge(data.table(EnsemblID = geneIDs), annot, by = "EnsemblID", all.x = T, sort = F)
            fit.eb$genes  <- annotation

            # 5.6) Output
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

write.csv(outputFileRegistry, paste0("RNASeq-pipeline_OutputRegistry.csv"), row.names = F)
