setwd("~/Sites/cats/sites/all/libraries/d3.fd/")
options(stringsAsFactors = FALSE)

real_fc <- function(x) {
  x=2^x
  for(i in 1:length(x))
    if(x[i]<1) {
      x[i]=-1/x[i] 
    }
  return(x)
}


# Read in data from gene mania
dat <- read.table("genemania-interactions-5.txt",  sep="\t", header=TRUE)
head(dat)

dat2=dat[dat$Network.group == "Physical Interactions" | dat$Network.group == "Pathway" |
         dat$Network.group == "Predicted" ,]

all=c(dat2$Gene.1,dat2$Gene.2)
id=unique(all)
group=rep(1,length(id))
dat1=cbind(id, group)

dat2$source = match(dat2$Gene.1, id)-1
dat2$target = match(dat2$Gene.2, id)-1


dat2=dat2[,c("Gene.1","Gene.2","source","target")]
colnames(dat2)<-c("Entity.1","Entity.2",	"source",	"target")
write.csv(dat2,"snp-links-cn-synapse.csv", row.names = F)

# Get the expression and GWiS data
exp=read.csv("ROSMAP_PFC_FPKM_CpDxLow_AD-NCI.csv", sep=",", header=T)
gwis=read.csv("igap_gwis.GWiS.Summary.csv", sep=",", header=T)
dat1=data.frame(dat1)

dat1$fc=exp[match(dat1$id,exp$GeneSymbol),"logFC"]
dat1$fc[is.na(dat1$fc)] <- 0
dat1$fc=real_fc(dat1$fc)

dat1$log10_exp=exp[match(dat1$id,exp$GeneSymbol), "adj.P.Val"]
dat1$log10_exp[is.na(dat1$log10_exp)] <- 0
dat1$log10_exp=-log10(dat1$log10_exp)
dat1$log10_exp[is.infinite(dat1$log10_exp)]<-0

dat1$log10_IGAP=gwis[match(dat1$id,gwis$Name),"Pval"]
dat1$log10_IGAP[is.na(dat1$log10_IGAP)] <- 0
dat1$log10_IGAP=-log10(dat1$log10_IGAP)
dat1$log10_IGAP[is.infinite(dat1$log10_IGAP)]<-0

dat1$group[dat1$log10_exp > 1.3]=1
dat1$group[dat1$log10_exp < 1.3]=2

write.csv(dat1,"snp-genes-cn-synapse.csv", row.names = F)


##old method
#exp=read.csv("/Users/sdas/Dropbox (Partners HealthCare)/MSBB/mongo/Yann/mongo_ROSMAP_PFC_CpDxLow_AD-NCI.csv")
# exp=read.csv("~/Dropbox (Partners HealthCare)/MSBB/Pipelines/Results/RNASeq/ROSMAP_PFC_FPKM_CpDxLow_AD-NCI.csv", sep="\t", header=T)
# all = exp[match(id,exp$GeneSymbol), c("ROSMAP_PFC_FPKM_CpDxAll_AD_PC_logFC","ROSMAP_PFC_FPKM_CpDxAll_AD_PC_adj.P.Val","IGAP_1_Pvalue","MayoEGWAS_eQTL_TCX_AD_Pvalue")]
# colnames(all) = c("logFC","adj.P.Val","IGAP_1_Pvalue","MayoEGWAS_eQTL_TCX_AD_Pvalue")
# all[is.na(all)] <- 0
# all$fc= real_fc(all$logFC)
# all$log10_exp=-log10(all$adj.P.Val)
# all$log10_IGAP=-log10(all$IGAP_1_Pvalue)
# all$log10_eQTL=-log10(all$MayoEGWAS_eQTL_TCX_AD_Pvalue)
# all$log10_exp[is.infinite(all$log10_exp)]<-0
# all$log10_eQTL[is.infinite(all$log10_eQTL)]<-0
# all$log10_IGAP[is.infinite(all$log10_IGAP)]<-0
# 
# genes=read.csv("snp-genes-microglia-base.csv")
# dat3=cbind(genes,all)
# dat3$group[dat3$log10_exp > 1.3]=1
# dat3$group[dat3$log10_exp < 1.3]=2
# 
# dat3=dat3[,-c(3,4,5,6)]
# write.csv(dat3,"snp-genes-microglia.csv", row.names = F)
# 


source("https://bioconductor.org/biocLite.R")
biocLite("STRINGdb")
library(STRINGdb)
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
# map genes to stringdb ids and graph genes
example1_mapped <- string_db$map( dat1, "id", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id
string_db$plot_network( hits )

ints=string_db$get_interactions( hits )

## Calcineurin genes
cn=read.table("~/Dropbox (Partners HealthCare)/DataLENS_Paper/Networks/Calcineurin/CN-genes.txt",  sep="\t", header=TRUE)
cn=unique(cn$Human_Gene)
neuron=read.table("~/Dropbox (Partners HealthCare)/DataLENS_Paper/Networks/Calcineurin/geneset.neuron.txt",  sep="\t", header=TRUE)
synapse=read.table("~/Sites/cats/sites/all/libraries/d3.fd/geneset.synapse.txt", sep="\t", header=TRUE)
nfat=read.table("~/Dropbox (Partners HealthCare)/DataLENS_Paper/Networks/Calcineurin/geneset.NFAT.txt",  sep="\t", header=TRUE)
creb=read.table("~/Sites/cats/sites/all/libraries/d3.fd/geneset.creb.txt",  sep="\t", header=TRUE)
mir24=read.table("~/Sites/cats/sites/all/libraries/d3.fd/geneset.mir24.txt",  sep="\t", header=TRUE)
barres=read.table("~/Dropbox (Partners HealthCare)/DataLENS_Paper/NFAT/barreslab_rnaseq.csv", sep=",", header=T)
nfat_up=read.table("~/Dropbox (Partners HealthCare)/DataLENS_Paper/NFAT/ROSMAP-UP-NFAT-Targets.txt",sep="\t", header=TRUE)
  
intersect(neuron$GO_NEURON_PART,synapse$GO_SYNAPSE)
cn_neuron=intersect(cn,neuron$GO_NEURON_PART)
write.csv(cn_neuron,"cn_neuron.csv", row.names = F)
intersect(cn_neuron,gwis$Name)

cn_nfat=intersect(cn,nfat$TGGAAA_NFAT_Q4_01)
write.csv(cn_nfat,"cn_nfat.csv", row.names = F)
intersect(cn_nfat,gwis$Name)

cn_synapse=intersect(cn,synapse$GO_SYNAPSE)
write.csv(cn_synapse,"cn_synapse.csv", row.names = F)
intersect(cn_synapse,gwis$Name)

intersect(neuron$GO_NEURON_PART,synapse$GO_SYNAPSE)
intersect(cn_neuron,synapse$GO_SYNAPSE)
cn_snypase_nfat=intersect(cn_synapse,nfat$TGGAAA_NFAT_Q4_01)
cn_synapse_creb=intersect(cn_synapse,creb$CREB_Q3)
intersect(cn,creb$CREB_Q3)
intersect(cn_synapse, mir24$CTGAGCC_MIR24)
nfat_exp=exp[exp$GeneSymbol %in% cn_snypase_nfat,]
nfat_exp1=exp[exp$GeneSymbol %in% cn_nfat,]

## Barres cell types
barres$Gene.symbol = toupper(barres$Gene.symbol)
nfat_up_barres=barres[match(nfat_up$Genes,barres$Gene.symbol),]
nfat_up_barres=na.omit(nfat_up_barres)
rownames(nfat_up_barres)=nfat_up_barres$Gene.symbol
barres_cells=nfat_up_barres[,3:9]
xx = data.frame(unlist(apply(barres_cells,1,function(x) colnames(barres_cells)[which(x==max(x))])))

xx=xx[rownames(xx) %in% nfat_up_barres$Gene.symbol,]
table(xx)
setdiff(rownames(xx), nfat_up_barres$Gene.symbol)
                                         

# CREB test
1-phyper(6, 140, 20000, 253)
