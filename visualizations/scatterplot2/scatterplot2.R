setwd("~/Sites/cats/sites/all/libraries/d3.scattercsv2")
options(stringsAsFactors = FALSE)
library(data.table)


dat <- fread("results.tsv")
head(dat)

dat1 = dat[,c("GeneSymbol","ROSMAP_PFC_FPKM_CpDxLow_AD_NCI_logFC","MayoEGWAS_eQTL_TCX_AD_Pvalue","IGAP_1_Pvalue")]
colnames(dat1) = c("GeneSymbol","logFC","color","Pval")
dat1$Pval = -log10(dat1$Pval)
dat1 = dat1[!is.na(dat1$Pval)]

range(dat1$color, na.rm = T)
dat1[is.na(dat1),] = 1
sum(dat1$color < .000001)
dat1$color[dat1$color > .000001] = "Not Significant"
dat1$color[as.numeric(dat1$color) < .000001] = "eQTL pval < 1.0E-6"
dat1$color[dat1$Pval >= 6] = "GWAS pval < 1.0E-6"
table(dat1$color)
write.table(dat1, file="ROSMAP-logFC-GWAS-Pval.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
