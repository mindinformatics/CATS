# control: 3 CKp25_2wk
# case: 3 CKp25_6wk
library(DESeq2)

setwd("~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/MouseResults/")

if(dir.exists("Data")==FALSE){
  dir.create("Data")
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65159&format=file&file=GSE65159%5FrnaSeqCount%2Etxt%2Egz", destfile = "Data/rnaseq.txt")
}

gset = read.table("Data/rnaseq.txt")

CKp25_2wk = c(10,11,12)
CKp25_6wk = c(4,5,6)

gset.selection<-gset[,c(CKp25_2wk, CKp25_6wk)] 

condition <- factor(c("Control","Control", "Control","Case","Case","Case"))
condition <- relevel(condition, ref = "Control")

dds<- DESeqDataSetFromMatrix(countData = gset.selection, 
                             colData = DataFrame(condition),
                             design = ~condition)

# filter rows with low counts to avoid spurious large FC's
#dds <- dds[ rowSums(counts(dds)) > 20, ]

# Tsai Method of filtering every row element > 20 reads
index = as.vector(which(apply(X = counts(dds), FUN = function(x) all(x>20), MARGIN = 1)))
dds <- dds[index, ]

dds.deseq<-DESeq(dds)
dds.results<-results(dds.deseq)

n<-dim(dds.results)[1]

write.csv(x = dds.results[1:(n-4), -1], 
          file = "GSE65159_Tsai_CKp25_6wkvsCKp25_2wk.csv", row.names = T)