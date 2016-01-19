#' Transform intensity values in a gradient of color from green to red.
#' @author Simon J Pelletier
#' @param expr.matrix A vector of the names of all samples
#' @param results A list of matrices. Each matrix contain a
#' @param threshold The lowest value of correlation to make the link
#' considered significantly high
#' @param typeID The type of ID used, e.g.
#' @return A vector of all possible comparisons
#' @keywords cytoscape export
#' @seealso \code{\link{}}
#' @examples
#' expr.matrix
#' results
#' networkSelectedComparison = export2cytoscape(expr.matrix,results[[selectedComparison]],threshold,typeID)

#' @export
export2cytoscape = function(expr.matrix,results,threshold=0.9,typeID="ensembl_gene_id"){
  expr.matrix.selected = expr.matrix[match(as.character(results[,typeID][!is.na(results[,typeID])]),rownames(expr.matrix)),]

  expr.adj.selected = adjacency(matrix(as.numeric(t(expr.matrix.selected)),ncol=nrow(expr.matrix.selected),nrow=ncol(expr.matrix.selected)))
  rownames(expr.adj.selected) = rownames(expr.matrix.selected)
  colnames(expr.adj.selected) = rownames(expr.matrix.selected)
  network = exportNetworkToCytoscape(expr.adj.selected, threshold = threshold)
  network$nodeData = merge(network$nodeData,results, by.x = "nodeName", by.y=typeID,all.y=TRUE)
  network$edgeData = cbind(network$edgeData[,1:4],rep("adjacency(WGCNA)",nrow(network$edgeData)))
  colnames(network$edgeData)[ncol(network$edgeData)] = "type"

  numChildren = rep(0,nrow(network$nodeData))
  #pourrait peut-etre etre changé pour utiliser with()
  for(n in 1:nrow(network$nodeData)){
    if (!is.na(as.character(network$nodeData[,1][n]))){
      numChildren[n] = length(grep(as.character(network$nodeData[,1][n]),c(as.character(network$edgeData[,1]),as.character(network$edgeData[,2]))))
    } else {
      numChildren[n] = 0
    }
  }
  network$nodeData = cbind(symbol=network$nodeData,numChildren)
  #filter="symbol"
  #write.csv(network$nodeData[,-12],file = paste0("node_rnaseq_",name,".csv"))
  # write All result GO
  #htmlReport(go[[i]] , file= paste0("ontology_", colnames(contrast.matrix)[i]), ".html")
  return(network)
}
