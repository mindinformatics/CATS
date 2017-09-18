source('~/Desktop/CATS/MousePipeline/Code/Download.R')
source('~/Desktop/CATS/MousePipeline/Code/CaseControl.R')
dat = DownloadFirstElement("GSE9566")

#check for whether log tranform needs to done
summary(exprs(dat)[1,])
pData(phenoData(dat))$`data_processing`[1]


###################################################################################
# Identify Cases and Controls

titles = pData(phenoData(dat))$title

defFactors<-function(title){
  #filter into different categories p indicates the post-natal days
  #Astro's P7-P8
  index<-grep(x=title,pattern="Astros P[1-8] ")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors<-data.frame(factor(tempvec))
  names(myFactors)<-"AstroP7"
 
   #Astro's P17
  index<-grep(x=title, pattern="Astros P[1][7] |Astros P[3][0]")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$AstroP17<-factor(tempvec)
  
  #Astro's P17-grey matter (P17g)
  index<-grep(x=title,pattern="Astros P17g ")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$AstroP17g<-factor(tempvec)
  
  #Neurons P7
  index<-grep(x=title,pattern="Neurons P[0-7] ")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$NeuronsP7 <-factor(tempvec)
  
  # Neurons P17
  index<-grep(x=title,pattern="Neurons P[1-3][0-9] ")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$NeuronsP17<-factor(tempvec)
  
  # Neurons endothelial cell depleted P7n P17n
  index<-grep(x=title,pattern="Neurons P[0-9][0-9]n ")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$NeuronsP7n<-factor(tempvec)
  
  # OPC's
  index<-grep(x=title,pattern="OPCs")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$OPCs<-factor(tempvec)
  
  # GalC-OL's
  #  galactocerebroside (GalC) first Oligodendrocyte specific marker to be expressed
  index<-grep(x=title,pattern="OLs.*GC_[a-z]")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$GalCOLs<-factor(tempvec)
  
  #MOG-OL's
  # MOG is a marker for mature Oligodendrocytes
  index<-grep(x=title,pattern="Myelin OLs.*MOG_[a-z]")
  tempvec<-rep(0, length(title))
  tempvec[index]<-1
  myFactors$MOGOLs<-factor(tempvec)
  
  return(myFactors)
}

Exprs1.factors<-defFactors(as.vector(titles))

Astros<-ifelse(Exprs1.factors$AstroP7==1|Exprs1.factors$AstroP17==1|Exprs1.factors$AstroP17g==1, 1, 0)
Neurons<-ifelse(Exprs1.factors$NeuronsP7==1|Exprs1.factors$NeuronsP17==1|Exprs1.factors$NeuronsP7n==1, 1, 0)
OLs<-ifelse(Exprs1.factors$OPCs==1|Exprs1.factors$MOGOLs==1|Exprs1.factors$GalCOLs==1, 1, 0)

CClist = CaseControl(Astros, dat)

