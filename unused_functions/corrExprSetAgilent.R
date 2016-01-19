#Not used at the moment
corrExprSetAgilent = function(expr.matrix) {
  expr.matrix = expr.matrix[!is.na(rownames(expr.matrix)),]
  expr.matrix = expr.matrix[grep("^A_",rownames(expr.matrix)),]
  for(i in 1:nrow(expr.matrix)){
    expr.matrix[i,] = apply(matrix(expr.matrix[grep(rownames(expr.matrix)[i],rownames(expr.matrix)),],ncol = ncol(expr.matrix)),2,mean)
  }
  return(expr.matrix)
}
