#THIS FUNCTION DOESN'T SEEM TO BE USED
changeID <- function(expr.matrix,typeID,genes_annotation_unique){
  symbols = as.character(genes_annotation_unique[match(rownames(expr.matrix),genes_annotation_unique[,"ID"]),typeID])
  pos <- !is.na(symbols) & !duplicated(symbols)
  expr.matrix2 <- expr.matrix[pos,]
  rownames(expr.matrix2) = symbols[pos]
  return(expr.matrix2)
}
