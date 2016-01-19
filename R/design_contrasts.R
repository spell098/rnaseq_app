#' Create a design matrix for linear models analysis
#' @author Simon J Pelletier
#' @import limma
#' @param names A vector of the names of all samples
#' @return a design matrix for linear models analysis
#' @keywords design linear limma
#' @seealso
#' \code{\link{design_contrasts}}
#' \code{\link{makeContrasts}}
#' \code{\link[limma]{lmFit}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' design <- design_contrasts(expr.matrix)
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' @export
design_contrasts = function(names){
  #names=strsplit(colnames(expr.matrix),"_")

  #names2=unlist(lapply(names,function(x){
  #  paste(x[-length(x)],collapse="_")
  #}))
  #colnames(expr.matrix) = names
  f <- factor(names, levels = unique(names))
  design <- model.matrix(~0+f)
  colnames(design) <- levels(f)
  rownames(design) = names
  return(design)
}
