#' Find to which genes belong to the same modules
#' @import RCytoscape
#' @author Simon J Pelletier
#' @param expr.matrix  A matrix of expression values. Rows are genes, columns are samples
#' @param ncores Number of cores available for parallel programming (foreach function)
#' @param names.unique Names of all groups of sample
#' @return A list of WGCNA results
#' @keywords cytoscape export
#' @seealso \code{\link{RCytoscape}}
#' @examples
#' ONLY USE A SUBSET; CAN TAKE A VERY LONG TIME TO RUN
#' results
#' networkSelectedComparison = export2cytoscape(expr.matrix,results[[selectedComparison]],threshold,typeID)
#' @export
WGCNA_modules = function(expr.matrix,ncores=1,names.unique){
  print("Running WGCNA to find modules")
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads(nThreads=ncores)
  nSets = length(names.unique);
  setLabels = names.unique
  # Form multi-set expression data: columns starting from 9 contain actual expression data.
  multiExpr = vector(mode = "list", length = nSets)

  for(i in 1:length(names.unique)){
    multiExpr[[i]] = list(data = as.data.frame(t(expr.matrix[,grep(names.unique[i],colnames(expr.matrix))])))
  }
  exprSize = checkSets(multiExpr)
  powers = c(seq(4,10,by=1), seq(12,20, by=2));
  powerTables = vector(mode = "list", length = nSets);
  # Call the network topology analysis function for each set in turn
  saveRDS(powerTables,file="powerTables.rds")

  print("Looking for soft thresholds")
  powerTables = foreach(set=1:nSets) %dopar% {
    #print(paste0("Soft threshold #",set," / ",nSets))
    list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                  verbose = 2)[[2]])
  }
  collectGarbage()
  print("Looking for consensus modules")
  bnet = blockwiseConsensusModules(
    multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25, numericLabels = TRUE,
    minKMEtoStay = 0,
    saveTOMs = TRUE, verbose = 5,
    nThreads=ncores)
  save(bnet,file = "bnet.rdata")
  print("Modules informations saved as bnet.rdata")
  return(bnet)
}
