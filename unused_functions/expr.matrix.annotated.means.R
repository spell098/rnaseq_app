#' UNFINISHED: the means are hard coded; change to select columns automatically line 33
#' Create a design matrix for linear models analysis
#' @author Simon J Pelletier
#' @import limma
#' @param expr.matrix.annotated A matrix of values with annototations (WHICH???) for each element (rows)
#' @param results topTable of significant genes
#' @return a design matrix for linear models analysis
#' @keywords design linear limma
#' @seealso
#' \code{\link{design_contrasts}}
#' \code{\link{makeContrasts}}
#' \code{\link[limma]{lmFit}}
#' \code{\link[limma]{topTable}}

#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' design <- design_contrasts(expr.matrix)
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' @export

expr.matrix.annotated.means = function(results,expr.matrix.annotated){
  col_clu=redgreen(100)
  col_clu = col_clu[100:1]

  signifGenes = unique(c(results[[1]]$symbol,results[[2]]$symbol,results[[3]]$symbol))
  signifGenes = signifGenes[signifGenes != ""]
  expr.matrix.annotated.signif = expr.matrix.annotated[match(signifGenes,expr.matrix.annotated$symbol),]
  rownames(expr.matrix.annotated.signif) = expr.matrix.annotated.signif$symbol
  expr.matrix.signif = data.matrix(expr.matrix.annotated.signif[,3:13])
  expr.matrix.signif.means = cbind(apply(expr.matrix.signif[,1:4],1,mean),apply(expr.matrix.signif[,5:8],1,mean),apply(expr.matrix.signif[,9:11],1,mean))
  colnames(expr.matrix.signif.means) = names(results)
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  return(expr.matrix.signif.means)
}
