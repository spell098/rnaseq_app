numberNames = function(names){ 
  names=as.character(names)
  names2=NULL
  for (i in unique(names)){
    l=grep(i,names)
    k=0
    for(j in l){
      k=k+1
      names2[j] = paste(names[j],k,sep="_")
    }
  }
  return(names2)
}

