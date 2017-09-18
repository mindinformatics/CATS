CaseControl = function(ContrastVector, GEOobject){
  # eliminate samples marked as "X"
  sel <- which(ContrastVector != "X")
  ContrastVector <- ContrastVector[sel]
  gset <- exprs(GEOobject)[ ,sel]
  
  return(list(gset, ContrastVector))
  }


