signifRows = function(expr.matrix,resultsSummary){
  expr.matrix.signif = data.matrix(expr.matrix[match(rownames(resultsSummary),rownames(expr.matrix)),])
  colnames(expr.matrix.signif) = colnames(expr.matrix)
  rownames(expr.matrix.signif) = as.character(resultsSummary[,1])
  return(expr.matrix.signif)
}
