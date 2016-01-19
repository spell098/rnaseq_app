#' Transform intensity values in a gradient of color from green to red.
#' @author Simon J Pelletier
#NOT USED ANYMORE APPARENTLY

#' @import GEOquery
#' @param comparisonSelection Selected comparison of 2 groups
#' @return A vector of colors. The darkest green values are the lowest and the darkest red values are the highest.
#' @keywords comparisons
#' @seealso
#' \code{\link[GEOquery]{getGEO}}
#' \code{\link[Biobase]{ExpressionSet}}
#' @examples
#' gset <- getGEO('GSE61276', GSEMatrix =TRUE) #GSE61276 GSE12654
#' exprset <- gset[[1]]
#' comparisons <- comparisonsPheno(exprset)[[1]]
#' comparisonsTable <- comparisonsPheno(exprset)[[2]]
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
