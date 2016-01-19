#' Add color to each sample. Samples with the same color belong to the same group
#' @author Simon J Pelletier
#' @param names Names of each sample
#' @param names.unique Name of each group
#' @return A vector of strings corresponding to the color of all samples.
#' @export
color_groups = function(names,names.unique){
  colors = vector("character",length(names))
  for (i in 1:length(names.unique)){
    colors[(grep(names.unique[i],names))] = wes()[i]
  }
  return(colors)
}
