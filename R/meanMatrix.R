meanMatrix = function(expr.matrix.signif,names.unique){
  expr.matrix.signif.means = matrix(ncol=length(names.unique),nrow=nrow(expr.matrix.signif))
  
  for (i in 1:length(names.unique)){
    expr.matrix.signif.means[,i] = data.matrix(apply(expr.matrix.signif[,grep(names.unique[i],colnames(expr.matrix.signif))],1,mean))
  }
  colnames(expr.matrix.signif.means) = names.unique
  rownames(expr.matrix.signif.means) = rownames(expr.matrix.signif)
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  return(expr.matrix.signif.means)
}
