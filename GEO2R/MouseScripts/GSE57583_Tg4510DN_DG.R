# control: 6 samples DN dentate gyrus
# case: 6 samples Tg4510 dentate gyrus

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Sep 18 15:39:04 EDT 2017
setwd("~/Desktop/CATS/GEO2R/GSE57583_BMS/")
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE57583", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8759", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "XXXXXXXXXXX000000XXXXXXXXXXX111111XXXXXXXXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10e10)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="Results/GSE57583_Tg4510DN_DG.csv", row.names=F, sep=",")


# ################################################################
# #   Boxplot for selected GEO samples
# library(Biobase)
# library(GEOquery)
# 
# # load series and platform data from GEO
# 
# gset <- getGEO("GSE57583", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL8759", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# # group names for all samples in a series
# gsms <- "XXXXXXXXXXX000000XXXXXXXXXXX111111XXXXXXXXXXXXXXX"
# sml <- c()
# for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
# sml <- paste("G", sml, sep="")  set group names
# 
# # eliminate samples marked as "X"
# sel <- which(sml != "X")
# sml <- sml[sel]
# gset <- gset[ ,sel]
# 
# # order samples by group
# ex <- exprs(gset)[ , order(sml)]
# sml <- sml[order(sml)]
# fl <- as.factor(sml)
# labels <- c("Control","Case")
# 
# # set parameters and draw the plot
# palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
# dev.new(width=4+dim(gset)[[2]]/5, height=6)
# par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
# title <- paste ("GSE57583", '/', annotation(gset), " selected samples", sep ='')
# boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
# legend("topleft", labels, fill=palette(), bty="n")