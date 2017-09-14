
CaseControl = function(GEOobject){
  Comparisons<-pData(phenoData(GEOobject))$`genotype/variation`
 
  APLP2ko_WT = ifelse(Comparisons=="APLP2 knockout", 1, ifelse(Comparisons=="wild type", 0, NA))
  APPko_WT = ifelse(Comparisons=="APP knockout", 1, ifelse(Comparisons=="wild type", 0, NA))
  APPsaki_WT = ifelse(Comparisons=="APPsa knockin", 1, ifelse(Comparisons=="wild type", 0, NA))
  
  return(data.frame(cbind(APLP2ko_WT, APPko_WT, APPsaki_WT)))
}
