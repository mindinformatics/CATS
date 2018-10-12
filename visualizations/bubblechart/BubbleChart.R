setwd("~/Sites/cats/sites/all/libraries/d3.bubblechartcsv")
options(stringsAsFactors = FALSE)
library(data.table)
library(dplyr)
currdate=format(Sys.Date(), format="%Y-%m-%d")

realfc <- function (x) {
  fc=ifelse(x>0,2^x,-1/2^x)
  return(fc)
}

# Directions:
# 1. Put gene symbols into input file: input_human_genes.txt
# 2. Select Contrast/Data Types below
# 3. Set output file name as data backup (just below selecting Contrast)


# Must use human gene symbols, one per line, header is GeneSymbol.
genes_list = read.csv("input_human_genes.txt", sep="\t", header = T)
genes = genes_list[,"GeneSymbol"]

analysisFile = read.table("~/Dropbox (Partners HealthCare)/MSBB/Methods/Output/main.csv", sep=",", header=T, stringsAsFactors=F)
head(analysisFile)
colnames(analysisFile)

# Select Contrast/StratFactor/DataType, etc
analysisFile=analysisFile[(analysisFile$StratFactor == "CpDxLow" & analysisFile$Contrast == "AD-NCI" & analysisFile$DataType != "RC"), ]

# Set Output File Name based on Contrast
fname2 = "Paper-CpDxLow-AD-NCI-Bubblechart-ZhangGenes.csv"


# Create output dataframe.
dat <- data.frame(study=character(),bregion=character(),cbregion=character(),dtype=character(),contrast=character(),
                  type=character(), GeneSymbol=character(), logFC=double(), FC=double(),
                  P.Value=double(), adj.P.Val=double(), stringsAsFactors=FALSE)

sapply(seq_along(analysisFile$FileName), function(i) {if (analysisFile$PipelineScript[i] == "uArray_pipeline.R" | analysisFile$PipelineScript[i] == "RNASeq_pipeline.R") {
  bregion = analysisFile$BrainRegionFull[i]
  cbregion = analysisFile$ClusteredBrainRegion[i]
  study = analysisFile$StudyName[i]
  contrast = analysisFile$Contrast[i]
  type=analysisFile$DataType[i]
  if (analysisFile$PipelineScript[i] == "uArray_pipeline.R") { dtype <- "Microarray" }
  if (analysisFile$PipelineScript[i] == "RNASeq_pipeline.R") { dtype <- "RNA-Seq" }
  fname= paste("~/Dropbox (Partners HealthCare)/MSBB/Pipelines/",analysisFile$Folder[i], analysisFile$FileName[i], sep="")
  analysis<-fread(fname)

  x<-grep("Gene.symbol",colnames(analysis))
  if(length(x) == 0) {colnames(analysis)[x] <- "GeneSymbol"}

  v1<-filter(analysis, GeneSymbol %in% genes)

  #v2<- v1[as.numeric(v1$P.Value)< 0.05,]
  v2<-do.call(rbind, lapply(split(v1,as.factor(v1$GeneSymbol)), function(x) {return(x[which.min(x$P.Value),])}))

  v3<-cbind(study, bregion, cbregion,  dtype, contrast, type, v2$GeneSymbol,v2$logFC, realfc(v2$logFC), v2$P.Value, v2$adj.P.Val)
  
  # Add results back; Note use of <<- assign back to global variable rather than use local copy within function
  v4 = as.data.frame(v3)
  colnames(v4) =c("study", "bregion", "cbregion", "dtype", "contrast", "type", "GeneSymbol", "logFC", "FC", "P.Value", "adj.P.Val")
  dat <<- rbind(dat,v4)
  check = paste0(dtype, ": ", study,": ",bregion,": ","Genes: ", dim(v4)[1])
  #print(check)
  return(check)
  }
}
)
# Error Checking
tot = dim(analysisFile)[1] * dim(genes_list)[1]
if (tot == dim(dat)[1]) { print("All genes found in all studies!")} else {print("Gene(s) Not found in 1 or more studies!")}

# Convert columns that got changed to character back to numeric for processing
# Else creating p1: Error in dat_br$P.Value * N : non-numeric argument to binary operator
dat[,c(8:11)] <- sapply(dat[, c(8:11)], as.numeric)

# TODO
# If there are multiple probes per gene, pick the probe with the highest significance (lowest Pvalue)
xx=data.frame(table(dat$study, dat$bregion, dat$dtype, dat$GeneSymbol)) ## check for multiple probes
#xx[xx$Freq >1,]
if (dim(xx[xx$Freq >1,])[1] > 0) { print("Error! Remove multiple probes!")}


## 2 Step Pvalue Correction
# First adjust for RNA-Seq and microarray for the 5 regions that have both
genes=unique(dat$GeneSymbol)
p1 = rbind.data.frame()
for(gene in genes) {
  print(gene)
  dat_gene=dat[dat$GeneSymbol == gene,]

  bregions = unique(dat_gene$bregion)
  for( br in bregions) {
    print(br)
    dat_br=dat_gene[dat_gene$bregion == br,]
    dat_br=dat_br[order(dat_br$P.Value),]
    N=dim(dat_br)[1]
    rank=1:nrow(dat_br)
    p1=rbind(p1,cbind(dat_br,P1.BF=dat_br$P.Value*N,P1.BH=dat_br$P.Value*N/rank))
  }
}

# Next adjust for 19 brain regions
p2 = rbind.data.frame()
for(gene in genes) {
  print(gene)
  dat_gene=p1[p1$GeneSymbol == gene,]
  dat_gene=dat_gene[order(dat_gene$P.Value),]

  bregions = unique(dat_gene$bregion)
  N=dim(dat_gene)[1]
  rank=1:nrow(dat_gene)
  p.adj=p.adjust(dat_gene$P1.BH, method="BH")
  p2=rbind(p2,cbind(dat_gene,P2.BF=dat_gene$P1.BF*N,P2.BH=dat_gene$P1.BH*N/rank, p.adj))
}

p2$P1.BF = ifelse(p2$P1.BF > 1, 1,p2$P1.BF)
p2$P1.BH = ifelse(p2$P1.BH > 1, 1,p2$P1.BH)
p2$P2.BF = ifelse(p2$P2.BF > 1, 1,p2$P2.BF)
p2$P2.BH = ifelse(p2$P2.BH > 1, 1,p2$P2.BH)

sum(p2$P2.BH < 0.10)
sum(p2$p.adj < 0.25)

# Create an input file for bubble chart
dat1 = p2[p2$p.adj < 0.2, c("study","bregion","dtype","contrast","GeneSymbol","logFC","FC","P.Value","p.adj") ]
colnames(dat1) = c("Study","parent","DataType","Contrast","name","LogFC","size","PValue","AdjPValue")

dat1$size = abs(log10(dat1$AdjPValue))

xx=data.frame(table(dat1$name, dat1$parent)) ## check for multiple probes

brparents = read.csv("BrainRegionParents.csv", header = T)

uparents = unique(dat1$parent)
dat2 = data.frame("", brparents[match(uparents, brparents$BrainRegion), "Parent"], "", "", uparents, "", "5", "","")
colnames(dat2)=c("Study","parent","DataType","Contrast","name","LogFC","size","PValue","AdjPValue")

uparents = unique(dat2$parent)
dat3 = data.frame("", brparents[match(uparents, brparents$BrainRegion), "Parent"], "", "", uparents, "", "5", "","")
colnames(dat3)=c("Study","parent","DataType","Contrast","name","LogFC","size","PValue","AdjPValue")

uparents = unique(dat3$parent)
dat4 = data.frame("", brparents[match(uparents, brparents$BrainRegion), "Parent"], "", "", uparents, "", "5", "","")
colnames(dat4)=c("Study","parent","DataType","Contrast","name","LogFC","size","PValue","AdjPValue")

dat_all=rbind(dat4, dat3, dat2, dat1)
colnames(dat_all)=c("Study","parent","DataType","Contrast","name","LogFC","size","PValue","AdjPValue")

# Abbreviations so that figure looks nice.  "Occipital lobe" removed manually.
dat_all = apply(dat_all, 2, function(y) gsub("Occipital Visual Cortex", "OVC", y, ignore.case = FALSE))
dat_all = apply(dat_all, 2, function(y) gsub("Anterior Cingulate", "Ant Cingulate", y, ignore.case = FALSE))
dat_all = apply(dat_all, 2, function(y) gsub("Middle Temporal Gyrus", "Mid Temp Gyrus", y, ignore.case = FALSE))
dat_all = apply(dat_all, 2, function(y) gsub("Limbic system", "Limbic", y, ignore.case = FALSE))


write.table(dat_all, fname2, sep =",", col.names = T, row.names = F, na="")
write.table(dat_all, "js-bubblechart-input.csv", sep =",", col.names = T, row.names = F, na="")


