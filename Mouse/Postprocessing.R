# Postprocessing total script
setwd("~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/")
library(data.table)
library(dplyr)

# read in necesary files
mouseAnalysis = readxl::read_excel("MouseAnalyses.xlsx")
mouseDatasets = readxl::read_excel("MouseDatasets.xlsx")
mouseFull = full_join(mouseAnalysis, mouseDatasets, by="GEO")
rm(mouseAnalysis, mouseDatasets)

MouseToHuman = fread("MouseToHuman/MouseToHumanGenes.csv")
MouseToHuman = MouseToHuman[!is.na(MouseToHuman$MouseGene) & !is.na(MouseToHuman$HumanGene),]

# denormalize files
addColumns = function(mouseRow){
  # Mouse denormalized columns: Filename, Brain Region, Perturbation, Data Type, Lab/Group, Study
  resultDataFrame = fread( paste("MouseResults/", mouseRow[7], sep="") )
  
  M = nrow(resultDataFrame)
  resultDataFrame[,"Filename"] = rep(mouseRow[7], M)
  resultDataFrame[,"BrainRegion"] = rep(mouseRow[24], M)  
  resultDataFrame[,"Perturbation"] = rep(mouseRow[9], M)
  resultDataFrame[,"DataType"] = rep(mouseRow[22], M)
  resultDataFrame[,"Group"] = rep(mouseRow[21], M)
  resultDataFrame[,"Study"] = rep(mouseRow[2], M)
  
  return(resultDataFrame)
} 

dontaddColumns = function(mouseRow){
  # Mouse denormalized columns: Filename, Brain Region, Perturbation, Data Type, Lab/Group, Study
  resultDataFrame = fread( paste("MouseResults/", mouseRow[7], sep="") )
  M = nrow(resultDataFrame)
  resultDataFrame[,"Filename"] = rep(mouseRow[7], M)
  return(resultDataFrame)
} 

# mouse genes to human genes
addFromMouseGene = function(resultFile){
  #fileName = paste("HumanGeneAdded/", resultFile$Filename[1], sep="")
  fileName = paste("Simple/", resultFile$Filename[1], sep="")
  resultFile$Gene.symbol = gsub("///.*$", "", x = resultFile$Gene.symbol)
  Index = match(resultFile$Gene.symbol, MouseToHuman$MouseGene)
  resultFile[,"HumanGene"] = MouseToHuman$HumanGene[Index]
  
  names(resultFile)[which(names(resultFile)=="ID")]= "ProbeID"
  names(resultFile)[which(names(resultFile)=="adj.P.Val")] = "adjPValue"
  names(resultFile)[which(names(resultFile)=="P.Value")] = "PValue"
  names(resultFile)[which(names(resultFile)=="Gene.symbol")] = "MouseGene"
  names(resultFile)[which(names(resultFile)=="Gene.title")] = "GeneTitle"
 
  fwrite(resultFile[!(resultFile$MouseGene==""), !"Filename"], file = fileName, row.names = F, 
         col.names = T, sep = ",")
}

# ensembl to human genes
addFromENS = function(resultFile){
  #fileName = paste("HumanGeneAdded/", resultFile$Filename[1], sep="")
  fileName = paste("Simple/", resultFile$Filename[1], sep="")
  resultFile$V1 = gsub("///.*$", "", x = resultFile$V1)
  Index = match(resultFile$V1, MouseToHuman$MUSENSG)
  resultFile[,"HumanGene"] = MouseToHuman$HumanGene[Index]
  resultFile[,"MouseGene"] = MouseToHuman$MouseGene[Index]
  
  names(resultFile)[which(names(resultFile)=="log2FoldChange")] = "logFC"
  names(resultFile)[which(names(resultFile)=="pvalue")] = "PValue"
  names(resultFile)[which(names(resultFile)=="padj")] = "adjPValue"
  names(resultFile)[which(names(resultFile)=="V1")] = "Ensembl"
  
  fwrite(resultFile[!(resultFile$MouseGene==""),], file = fileName, 
         row.names = F, col.names = T, sep = ",")
}

#denormalized = apply(X = mouseFull, MARGIN = 1, FUN = addColumns)
notdenormalized = apply(X = mouseFull, MARGIN = 1, FUN = dontaddColumns)

lapply(X=notdenormalized[c(-30:-32)], FUN = addFromMouseGene)
lapply(X=notdenormalized[c(30:32)], FUN = addFromENS)

#lapply(X=denormalized[c(-30:-32)], FUN = addFromMouseGene)
#lapply(X=denormalized[c(30:32)], FUN = addFromENS)
