cpg.annotation = function(cpg,annotation){
  list_cpg=vector("list",nrow(cpg))
  for (i in 1:nrow(cpg)){
    for (j in 1:nrow(annotation)){
      if ( (cpg[i,2]<annotation[j,3] && cpg[i,3]>annotation[j,2]) || (cpg[i,2]>annotation[j,3] && cpg[i,3]<annotation[j,2]) && as.character(cpg[i,1]) == as.character(annotation[j,1])){
        test=1
        if(is.null(list_cpg[[i]])){
          n=1
        }
        list_cpg[[i]][n] = as.character(annotation[j,4])
        n=n+1
      }
    }
  }
  return(list_cpg)
}
