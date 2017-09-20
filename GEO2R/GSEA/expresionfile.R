MakeGSEAtxt = function(gset, name){
  vals = data.frame(exprs(gset))
  vals[,'Description'] = as.vector(rep(NA,dim(vals)[1]))
  write.table( x = vals[,c(dim(vals)[2], 1:dim(vals)[2]-1)], 
        file = paste("GSEA/",name), row.names = T, sep = "\t")
}

MakeGSEAcls = function(gsms, name){
  
  write(gsms, file = paste(name,".txt"))
}