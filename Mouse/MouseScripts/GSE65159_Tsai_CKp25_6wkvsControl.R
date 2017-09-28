# control: 3 CK_6wk
# case: 3 CKp25_6wk
library(DESeq2)

setwd("~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/MouseResults/")

if(dir.exists("Data")==FALSE){
  dir.create("Data")
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65159&format=file&file=GSE65159%5FrnaSeqCount%2Etxt%2Egz", destfile = "Data/rnaseq.txt")
}

gset = read.table("Data/rnaseq.txt")

Control = c(1,2,3,7,8,9)
CKp25_6wk = c(4,5,6)

gset.selection<-gset[,c(Control, CKp25_6wk)] 

#need to set the rownames inorder to identify

condition <- factor(c("Control","Control", "Control","Control",
                      "Control", "Control","Case","Case","Case"))

condition <- relevel(condition, ref = "Control")

dds<- DESeqDataSetFromMatrix(countData = gset.selection, 
                             colData = DataFrame(condition),
                             design = ~condition)

# filter rows with low counts to avoid spurious large FC's
#dds <- dds[ rowSums(counts(dds)) > 20, ]

#Method of Filtering as Described in Paper
index = as.vector(which(apply(X = counts(dds), FUN = function(x) all(x>20), MARGIN = 1)))
dds <- dds[index, ]

dds.deseq<-DESeq(dds)
dds.results<-results(dds.deseq)

n<-dim(dds.results)[1]

write.csv(x = dds.results[1:(n-4), -1], 
          file = "GSE65159_Tsai_CKp25_6wkvsControl.csv", row.names = T)