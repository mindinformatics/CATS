setwd("~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/")
library(data.table)
library(dyplr)

setwd("Simple/")
MouseResultFiles = list.files(".",pattern = "*.csv")

MouseToHuman = fread("../Modules/MouseToHuman/MouseToHumanGenesNew.csv", header =T)

limitBest = function(mouseResultFile, alpha=0.05, top=100) {
   tempFile = fread(mouseResultFile, header = T, sep =",")
   tempFile[,"abslogFC"] = abs(tempFile$logFC)
   tempFile %>% filter(adjPValue<alpha) %>% 
     top_n(n=top, wt=abslogFC) %>%
     arrange(desc(abslogFC)) -> tempFile.sig
 
   # Add mouse Ensembl and human ensembl
   index = match(tempFile$MouseGene, MouseToHuman$geneMouse)
   
   tempFile[, "EnsemblHuman"] = MouseToHuman$EnsemblHuman[index]
   tempFile[, "EnsemblMouse"] = MouseToHuman$EnsemblMouse[index]
   
   newName = paste( gsub(".csv", "", mouseResultFile), "_top100.csv", sep="")
   sigFileName = paste("../Top100adjPVal0.05/", newName, sep="")
   
   fwrite(x = tempFile.sig[,-length(tempFile.sig)], file = sigFileName)
}

lapply(X=MouseResultFiles, FUN = function(x) limitBest(x))


#check output

setwd("../Top100adjPVal0.05/")
OutputFiles = list.files(".", pattern = "*.csv")

for (fileName in OutputFiles) {
 TempFile = fread(file = fileName )
 if (max(TempFile$adjPValue)>0.05){
   stop("Error: max adjPval>0.05")
 }
  
}
