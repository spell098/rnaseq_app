#' Mean matrix
#' @description For every row of the initial expression matrix, the mean for each group is calculated.
#' @return Matrix
#' \describe{
#'  \item{columns}{Groups (mean values)}
#'  \item{rows}{Genes}
#' }
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#' topTable3 <- readRDS("data/results_LGVD.rds")
#' resultsSummary <- results_summary(results,topTable3,adjust="BH",typeID="ensembl_gene_id")
#' expr.matrix.signif <- signifRows(expr.matrix,resultsSummary)
#'
#' @seealso
#' \code{\link{signifRows}}
#' @export
meanMatrix = function(expr.matrix.signif){
  names.unique = unique(get_names(expr.matrix.signif))
  expr.matrix.signif.means = matrix(ncol=length(names.unique),nrow=nrow(expr.matrix.signif))

  for (i in 1:length(names.unique)){
    expr.matrix.signif.means[,i] = data.matrix(apply(expr.matrix.signif[,grep(names.unique[i],colnames(expr.matrix.signif))],1,mean))
  }
  colnames(expr.matrix.signif.means) = names.unique
  rownames(expr.matrix.signif.means) = rownames(expr.matrix.signif)
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  return(expr.matrix.signif.means)
}
