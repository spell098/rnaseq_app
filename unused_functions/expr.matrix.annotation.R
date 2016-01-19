#' Merge two matrices that have an agilentID column
#' @author Simon J Pelletier
#' @param ID.annotated.unique Matrix with a variety of annotations
#' @param expr.toBind Expression matrix annotated with modules numbers
#' @return
#' \describe{
#'  \item{expr.matrix.annotated}{Matrix annotated with all IDs (WHICH???)}
#' }
#' @keywords expr.toBind
#' @examples
#'
#' @export
expr.matrix.annotation = function(ID.annotated.unique,expr.toBind){
  colnames(expr.toBind)[1:2] = c("module","agilentID")
  expr.matrix.annotated = merge(expr.toBind,ID.annotated.unique,by.x = "agilentID", by.y = "efg_agilent_wholegenome_4x44k_v1")
  return(expr.matrix.annotated)
}
