#' Select in an expression matrix only the rows with significant elements (genes,transcripts...)
#' @author Simon J Pelletier
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param resultsSummary A matrix of data from every gene that is significanlty different
#' @return The significant rows of the initial expression matrix
#' between at least two groups. Columns = samples, rows = genes,transcripts,CpG...
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#' topTable3 <- readRDS("data/topTable3_LGVD.rds")
#' resultsSummary <- results_summary(results,topTable3)
#' @seealso
#' \code{\link[limma]{topTable}}
#' \code{\link{results_summary}}
#' @export
signifRows = function(expr.matrix,resultsSummary){
  expr.matrix.signif = data.matrix(expr.matrix[match(rownames(resultsSummary),rownames(expr.matrix)),])
  colnames(expr.matrix.signif) = colnames(expr.matrix)
  rownames(expr.matrix.signif) = as.character(resultsSummary[,1])
  return(expr.matrix.signif)
}
