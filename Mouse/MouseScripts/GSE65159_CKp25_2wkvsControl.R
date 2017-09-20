# control: 3 CK_2wk
# case: 3 CKp25_2wk
library(DESeq2)

if(dir.exists("Data")==FALSE){
  dir.create("Data")
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE65159&format=file&file=GSE65159%5FrnaSeqCount%2Etxt%2Egz", destfile = "Data/rnaseq.txt")
}

gset = read.table("Data/rnaseq.txt")

Control = c(1,2,3,7,8,9)
CKp25_2wk=c(10,11,12)

gset.selection<-gset[,c(Control, CKp25_2wk)] 

condition <- factor(c("Control","Control", "Control","Control",
                      "Control", "Control","Case","Case","Case"))
condition <- relevel(condition, ref = "Control")

dds<- DESeqDataSetFromMatrix(countData = gset.selection, 
                             colData = DataFrame(condition),
                             design = ~condition)

# filter rows with low counts to avoid spurious large FC's
dds <- dds[ rowSums(counts(dds)) > 20, ]
dds.deseq<-DESeq(dds)
dds.results<-results(dds.deseq)

n<-dim(dds.results)[1]

write.csv(x = dds.results[1:(n-4), -1],
          file = "Results/GSE65159_CKp25_2wkvsControl.csv")
