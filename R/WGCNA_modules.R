#' Find to which genes belong to the same modules
#' @import WGCNA
#' @author Simon J Pelletier
#' @param expr.matrix  A matrix of expression values. Rows are genes, columns are samples
#' @param ncores Number of cores available for parallel programming (foreach function)
#' @param names.unique Names of all groups of sample
#'
#' @return A list of WGCNA results
#' @keywords cytoscape export
#' @seealso \code{\link{WGCNA}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' results <- readRDS("data/results_LGVD.rds")
#'
#' bnet = WGCNA_modules(expr.matrix.sample)
#' @export
WGCNA_modules = function(expr.matrix){
  names.unique=unique(get_names(expr.matrix))
  print("Running WGCNA to find modules")
  nSets = length(names.unique)
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

  for(set in 1:nSets){
    print(paste0("Soft threshold #",set," / ",nSets))
    powerTables[[set]] = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                 verbose = 2)[[2]]

  }
  saveRDS(powerTables,"powerTables.rds")
  collectGarbage()
  print("Looking for consensus modules")
  bnet = blockwiseConsensusModules(
    multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25, numericLabels = TRUE,
    minKMEtoStay = 0,
    saveTOMs = TRUE, verbose = 5,
    nThreads=ncores)
  print("Modules informations saved as bnet.rdata")
  return(bnet)
}
