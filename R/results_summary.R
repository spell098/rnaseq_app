#' Select in an expression matrix only the rows with significant elements (genes,transcripts...)
#' @author Simon J Pelletier
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param resultsSummary A matrix of data from every gene that is significanlty different
#' @return The significant rows of the initial expression matrix
#' between at least two groups. Columns = samples, rows = genes,transcripts,CpG...
#' @examples
#' results <- readRDS("data/results_LGVD.rds")
#' topTable3 <- readRDS("data/topTable3_LGVD.rds")
#' resultsSummary <- results_summary(results,topTable3)
#' @seealso
#' \code{\link[limma]{topTable}}
#' @export
results_summary = function(results,topTable3,adjust="no",typeID="ensembl_gene_id"){
  genes = NULL
  if(adjust=="no") sortedBy = "P.Value" else sortedBy = "adj.P.Val"
  for(i in 1:length(results)){
    if(!is.null(results[[i]])){
      genes = unique(c(genes,as.character(results[[i]][,"symbol"])))
    }
  }
  print(genes)
  if(length(genes) > 0){
    if(typeID != "symbol"){
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,c("symbol","ensembl_gene_id","entrezgene_id")]
      rownames(resultsSummary) = resultsSummary[,typeID]
    } else {
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,c("symbol","entrezgene_id","ID")]
      rownames(resultsSummary) = resultsSummary[,typeID]
    }
    resultsSummary = resultsSummary[,-2]
    for (i in 1:length(topTable3)){
      resultsSummary[,i+2] = topTable3[[i]][match(genes,topTable3[[i]]$symbol),][,sortedBy]
      colnames(resultsSummary)[i+2] = names(topTable3)[i]
    }
  } else {
    if(typeID != "symbol"){
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,c("symbol","ensembl_gene_id","entrezgene_id")]
    } else {
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,typeID]
    }
  }
  return(resultsSummary)
}
