mean_duplicated = function(){
  for (i in 1:length(expr.matrix.annotated$ensembl_transcript_id)){
    x = grep(expr.matrix.annotated$ensembl_transcript_id[i],expr.matrix.annotated$ensembl_transcript_id)
    selected = as.matrix(expr.matrix.annotated[x,][,1:ncol(expr$A)+1])
    selected = matrix(as.numeric(selected),ncol = ncol(expr$A))
    mean = apply(selected,2,mean)
    mean.rbind = rbind(mean,mean,mean)
    expr.matrix.annotated[x,1:ncol(expr$A)+1] = mean
  }
  return(expr$A)
}
