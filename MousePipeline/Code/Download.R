source("http://bioconductor.org/biocLite.R")
library(GEOquery, quietly = T)

DownloadFirstElement = function(GEOtag) {
  GEOobject<-getGEO(GEO=GEOtag, GSEMatrix=T, AnnotGPL=T,
                    destdir= "~/Desktop/CATS/MousePipeline/Download/")
  
  if(length(GEOobject)!=1){GEOobject = GEOobject[[1]]}
  AddDatasetEntry(GEOobject, GEOtag)
  return(GEOobject)
}


DownloadFunction = function(GEOtag){
  GEOobject<-getGEO(GEO=GEOtag, GSEMatrix=T, AnnotGPL=T,
                    destdir= "~/Desktop/MousePipeline/Download/")
                      #"~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/Datasets/")
 
  #Cant ignore multipule platforms for multipule  
  #if (length(GEOobject)!=1) {
  #GEOobject = GEOobject[[1]]
  #}

  lapply(X = GEOobject, FUN = function(x) AddDatasetEntry(x,GEOtag)) 
  AddDatasetEntry(GEOobject, GEOtag)
  return(GEOobject)
}

AddDatasetEntry = function(GEOobject, GEOtag) {
  # Name, StudyLab, DataType, GEO, SampleSize, BrainRegion, Peturbations, TimeCourse, DataLinks
  # Species, Cell Type, OtherTissues, PubMedID, Description
  
  StudyLab = expinfo(experimentData(GEOobject))[['lab']]
  DataType = "" # every GPL needs to map to microarray or RNAseq 
  SampleSize = length(sampleNames(phenoData(GEOobject)))
  BrainRegion = paste(levels(pData(phenoData(GEOobject))$source_name_ch1), collapse = "//")
  Peturbations = ""
  TimeCourse = ""
  DataLinks = ""
  Species = paste(levels(pData(phenoData(GEOobject))$organism_ch1), sep=" ", collapse = "//") 
  CellType = ""
  OtherTissues = ""
  PubMedID = pubMedIds(experimentData(GEOobject))
  Platform = annotation(GEOobject)
  Description = ""
  
  Name = paste(GEOtag, StudyLab, sep = " ")
  
  Append = c(Name, StudyLab, DataType,GEOtag, SampleSize, BrainRegion, 
             Peturbations, TimeCourse, DataLinks,
             Species, CellType, OtherTissues, PubMedID, Platform, Description)
 
   
  #"~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/Datasets/Datasets.csv",
  Dataset = read.csv("~/Desktop/MousePipeline/Download/Datasets.csv",
                     header =T)
  
  # if(GEOtag %in% Dataset$GEO){}
  # have multipule datasets per GEO, deserve their own ID?
  Dataset[nrow(Dataset)+1, ] = Append
  write.csv(Dataset,"~/Desktop/MousePipeline/Download/Datasets.csv",
            row.names = F)
}
