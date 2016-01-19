namesSelectedComparisons = function(namesTable){
  namesSelected <- apply(namesTable,1,function(x){
    y<-lapply(x,function(z){
      paste(strsplit(z," ")[[1]],collapse="_")
    })
    paste(y,collapse="_")
  })
  names(namesSelected) = rownames(namesTable)
  return(namesSelected)
}
