#' @export
selectTopReaction = function(selectedComparison,topReactions){
  x = strsplit(selectedComparison,"-")[[1]]
  alt1=paste(x[1],x[2],sep="-")
  alt2=paste(x[2],x[1],sep="-")
  if (!is.null(topReactions[[alt1]])){
    selectedTopReactions = topReactions[[alt1]]
  } else {
    selectedTopReactions = topReactions[[alt2]]
  }
  return(selectedTopReactions)
}
