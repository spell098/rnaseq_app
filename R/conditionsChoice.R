#' Transform intensity values in a gradient of color from green to red.
#' @author Simon J Pelletier
#' @import GEOquery
#' @param comparisonSelection Selected comparison of 2 groups
#' @return A vector of colors. The darkest green values are the lowest and the darkest red values are the highest.
#' @keywords comparisons
#' @seealso
#' \code{\link[GEOquery]{getGEO}}
#' \code{\link[Biobase]{ExpressionSet}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' comparisonSelection <- "DORSAL\\VENTRAL"
#' names <- get_names(expr.matrix)
#' conditionsChoice(comparisonSelection,names)
#' @export
conditionsChoice = function(comparisonSelection,names){
  if (length(comparisonSelection) > 0){
    names2=NULL
    for (m in length(comparisonSelection):1){
      conditionsSelected = strsplit(comparisonSelection[m],'\\\\')
      for (j in conditionsSelected[[1]]){
        range=grep(j,names)
        names2[range] = j
      }
    }
  }
  return(names2)
}
