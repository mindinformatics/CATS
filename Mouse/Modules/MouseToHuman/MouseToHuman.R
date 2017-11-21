setwd("~/Dropbox (Partners HealthCare)/MSBB/Mouse 2.0/MouseToHuman/")

source("http://bioconductor.org/biocLite.R")
library(data.table)
library(dplyr)
library(annotationTools)

# Download full HomoloGene data file from ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current
# columns 
  #1) HID, homologene group ID
  #2) Taxonomy ID
  #3) Gene ID
  #4) Gene Symbol
  #5) Protein gi
  #6) Protein accession
homologene = read.delim("homologene.data", header=FALSE)

# Get an Inner Join of homologene from Mouse to Human
# 10090  Mus musculus
# 9606 Human
homologene.hs = homologene[homologene$V2==9606,]
names(homologene.hs) = c("HID", "species", "ncbiHuman", "geneHuman", "proteinGIHuman", "accessionHuman")
homologene.mm = homologene[homologene$V2==10090,]
names(homologene.mm) = c("HID", "species", "ncbiMouse", "geneMouse", "proteinGIMouse", "accessionMouse")

MouseToHuman = inner_join(homologene.hs, homologene.mm, by="HID")

# read in Human Gene information which includes the Human Symbol with an assocaited ncbiID
HumanGene = fread("ncbiGENE11.21.17/Homo_sapiens.gene_info")
MouseGene = fread("ncbiGENE11.21.17/Mus_musculus.gene_info")

indexHuman = match(MouseToHuman$ncbiHuman, HumanGene$GeneID)
indexMouse = match(MouseToHuman$ncbiMouse, MouseGene$GeneID)

HumanGene.matched = HumanGene[indexHuman, ]
MouseGene.matched = MouseGene[indexMouse, ]

MouseToHuman[,"EnsemblHuman"] = stringr::str_extract(string = HumanGene.matched$dbXrefs, "ENSG[0-9]+")
MouseToHuman[,"EnsemblMouse"] = stringr::str_extract(string = MouseGene.matched$dbXrefs, "ENSMUSG[0-9]+")

chosen = c("ncbiHuman", "geneHuman", "proteinGIHuman", "accessionHuman", "ncbiMouse",
           "geneMouse", "proteinGIMouse","accessionMouse","EnsemblHuman", "EnsemblMouse")
write.csv(MouseToHuman[,chosen], file = "MouseToHumanGenes.csv", row.names = F)
