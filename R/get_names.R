#' Get names of each sample and group names
#' @author Simon J Pelletier
#' @param names Samples names
#' @return return two objects:names1 -> all names without the number
#' (if their was one) which makes the name of samples in the same groups identical.
#' names.unique -> A vector of the name of each group (remove identical names from names1)
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' pval <- 0.05
#' summaryClust <- clusteringSummary(expr.matrix,pval)
#' @keywords names names.unique
#' @export
get_names = function(names){
  names1 = rep("",length(names))
  for (i in 1:length(names)){
    nameSplit = strsplit(names[i],"_")
    for (j in 1:(length(nameSplit[[1]]))){
      if(j==1) names1[i] = nameSplit[[1]][j]
      else if (is.na(as.numeric(gsub("([0-9]+).*$", "\\1", nameSplit[[1]][j])))) names1[i] = paste0(names1[i],"_",nameSplit[[1]][j])
    }
  }
  names.unique = unique(names1)
  return(list(names1,names.unique))
}
