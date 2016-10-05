# Last version of this script written on 25/SEP/2016
# by Yann Lambert, MIND Informatics
# ------------------

setwd("Desktop/AMP-AD/")

library(data.table)
library(dplyr)

rm(list = ls())

d <- fread("Report/Methods/Datasets.csv")
rf <- fread("Report/Methods/resultFiles.csv")

# For all result files get the corresponding dataset information
main <- data.table(merge(rf, d, by = "DatasetCode", all = T))

# Add for the GEO files the corresponding scripts
GEOScripts <- readLines("Report/Methods/GEOScripts.txt")
main[SourceType == "GEO"]$ProcessScript <- GEOScripts
main[SourceType == "GEO"]$PipelineScript <- GEOScripts


# ScriptDesc
scriptDesc <- list(
    Limma = "Analysis was performed by <a href=\"http://bioconductor.org/packages/release/bioc/html/limma.html\">Limma package</a>. No prior filtering was performed, neither on variance nor intensity.",
    RC = "TMM pre-normalization of Raw Counts was performed with  <a href=\"https://bioconductor.org/packages/release/bioc/html/edgeR.html\">EdgeR package</a>.",
    CPM = "Counts per million (CPM) were directly analyzed without further pre-normalization.",
    FPKM = "FPKM counts were directly analyzed with Limma with a log(x + 0.05) transformation."
)

uArrayAnnotDesc <- list(
    HBTRC = "The probe annotation file was downloaded from the GEO platform webpage (see above for reference)",
    MSBB = "The probe annotation file was downloaded from the GEO platform webpage (see above for reference)",
    ROSMAP = "<a href=\"https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv3.db.html\">illuminaHumanv3.db</a> package was used to select and remove poor quality Illumina mRNA probes. The same package was used for annotation.",
    MayoEGWAS = "<a href=\"https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanWGDASLv3.db.html\">illuminaHumanWGDASLv3.db</a> package was used to select and remove poor quality Illumina mRNA probes. The same package was used for annotation."
)

RNASeqAnnotDesc <- "Genes were annotated with the latest data available for GRCh38.p7 from Ensembl Biomart (http://useast.ensembl.org/biomart/martview/)"

main$ScriptDesc <- sapply(seq_along(main$FileName), function(i) {
    extra = ""
    annot = ""
    if(main[i]$DataType == "uArray" & main[i]$SourceType == "Synapse") {
        annot <- uArrayAnnotDesc[[main[i]$StudyName]]
    }
    else if(main[i]$DataType == "RC") {
        extra <- scriptDesc$RC
        annot <- RNASeqAnnotDesc
    }
    else if(main[i]$DataType == "CPM"){
        extra <- scriptDesc$CPM
        annot <- RNASeqAnnotDesc
    }
    else if(main[i]$DataType == "FPKM"){
        extra <- scriptDesc$FPKM
        annot <- RNASeqAnnotDesc
    }
    
    return(paste(scriptDesc$Limma, extra, annot, sep = " "))
})

# StratDesc
stratDesc <- list(
    CDR = paste0(readLines("Report/Methods/TextFiles/CDR.txt"), collapse = ""),
    Braak = paste0(readLines("Report/Methods/TextFiles/Braak.txt"), collapse = ""),
    MMSE = paste0(readLines("Report/Methods/TextFiles/MMSE.txt"), collapse = ""),
    CERAD = paste0(readLines("Report/Methods/TextFiles/CERAD.txt"), collapse = ""),
    CpDxStrict = paste0(readLines("Report/Methods/TextFiles/CompositeDiagnosis.txt"), collapse = ""),
    CpDxLow = paste0(readLines("Report/Methods/TextFiles/CompositeDiagnosis.txt"), collapse = ""),
    CpDxAll = paste0(readLines("Report/Methods/TextFiles/CompositeDiagnosis.txt"), collapse = ""),
    ClinicalDiagnosis = paste0(readLines("Report/Methods/TextFiles/ClinicalDiagnosis.txt"), collapse = "")
)

main$StratDesc <- sapply(seq_along(main$FileName), function(i) {
    return(stratDesc[[main[i]$StratFactor]])
})


fields <- list (
    "BrainRegionFull",
    "BrainRegionCode",
    "BMArea",
    "ClusteredBrainRegion",
    "DataType",
    "PlatformName",
    "PlatformGPL",
    "DataFile",
    "DataFileSynID",
    "CovFile",
    "CovFileSynID",
    "KeyFile",
    "KeyFileSynID",
    "ProcessScript",
    "PipelineScript",
    "ScriptDesc",
    "StratFactor",
    "StratDesc",
    "Contrast",
    "NSubjects",
    "NExcludedSubjects",
    "SubjectDistribution",
    "RefPubmedID"
)

template <- paste0(readLines("Report/Methods/TextFiles/Template.txt"), collapse = "")

main$MethodsHTML <- sapply(seq_along(main$FileName), function(i) {
        text <- template
        lapply(fields, function(x){
            pattern <- paste0("::",x,"::")
            text <<- sub(pattern, paste0(main[i][[x]]), text)
        })
        return(text)
    }
)

write.csv(main, "main.csv", row.names = F)

# write(text, "test.html")