library(data.table)

GEO2R = fread("raw.csv", header = T)
MYoutput = fread("GSE9566_FB_BARRES_ARRAYEXP_ASTROS.csv", header = T) 

# all the probes are ordered in the output the same
table(GEO2R$ID == MYoutput$V1)

# check that all numerical differences are small (likely due to rounding)
hist(MYoutput$logFC-GEO2R$logFC)
hist(MYoutput$PValue-GEO2R$P.Value)
hist(MYoutput$t-GEO2R$t)
hist(MYoutput$adjPValue-GEO2R$adj.P.Val)

# check for outliers, if anything is reasonably different 
range(MYoutput$logFC-GEO2R$logFC)
range(MYoutput$PValue-GEO2R$P.Value)
range(MYoutput$t-GEO2R$t)
range(MYoutput$adjPValue-GEO2R$adj.P.Val)
