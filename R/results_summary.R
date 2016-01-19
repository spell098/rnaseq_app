results_summary = function(results,topTable3,adjust,typeID){
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
