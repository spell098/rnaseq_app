#' Calculate results
#' @author Simon J Pelletier
#' @import limma
#' @param lm2.contrast .
#' @param expr.toBind .
#' @param pvalue .
#' @param logFC .
#' @param typeID .
#' @param genes_annotation_unique .
#' @param adjust .
#' @param annotations .
#' @return Results .
#' @keywords limma linear
#' @seealso
#' \code{\link[limma]{topTable}}
#' \code{\link{transcript_count}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' design <- design_contrasts(expr.matrix)
#' results_list <- results_topTable()
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' @export
results_topTable = function(lm2.contrast,expr.toBind,pvalue,logFC,typeID,genes_annotation_unique,annotations,adjust){
  results = vector("list", ncol(lm2.contrast$coefficients))
  topTable2 = vector("list", ncol(lm2.contrast$coefficients))
  topTable3 = vector("list", ncol(lm2.contrast$coefficients))
  names(results) = colnames(lm2.contrast$coefficients)
  names(topTable2) = colnames(lm2.contrast$coefficients)
  names(topTable3) = colnames(lm2.contrast$coefficients)

  for (i in 1:ncol(lm2.contrast$coefficients)){
    toptable = topTable(lm2.contrast , coef=i ,  sort.by = "none",
                        n=nrow(lm2.contrast$coefficients),adjust.method=adjust) #ENSG00000124831
    #hist(topTable)
    if (!is.null(expr.toBind)) {
      if (colnames(toptable)[1] != "ID"){
        topTable2[[i]] = cbind(expr.toBind[,"module"],rownames(toptable),toptable)
      } else {
        topTable2[[i]] = cbind(expr.toBind[,"module"],toptable)
      }
      colnames(topTable2[[i]])[1:2] = c("module",typeID)
    } else {
      if (colnames(toptable)[1] != "ID"){
        topTable2[[i]] = cbind(rownames(toptable),toptable)
      } else {
        topTable2[[i]] = toptable
      }
      colnames(topTable2[[i]])[1] = typeID
    }
    #print(head(topTable2[[i]]))
    topTable3[[i]] = merge(genes_annotation_unique,topTable2[[i]],by=typeID)
    #go_results[[i]] = merge(genes_ontology,)
    if(adjust=="no") sortedBy = "P.Value" else sortedBy = "adj.P.Val"

    mask = (topTable3[[i]]$logFC > logFC[2] | topTable3[[i]]$logFC < logFC[1]) & topTable3[[i]][,sortedBy] < pvalue
    results[[i]] = topTable3[[i]][mask,]

    ratio_transcript_signif = vector("numeric",nrow(results[[i]]))
    transcript_total = vector("numeric",nrow(results[[i]]))

    #pas efficace, mais ca marche  # Il y a des replicats de transcripts, pas bon... à fixer
    if(nrow(results[[i]]) >= 1){
      count_transcript = transcript_count(results[[i]],annotations,typeID)
      results[[i]] = cbind(results[[i]],count_transcript)
    }
    x <- data.frame(transcript_signif = rep(0,nrow(topTable3[[i]])),ratio_transcript_signif = rep(0,nrow(topTable3[[i]])))
    topTable3[[i]] = cbind(topTable3[[i]],x)

    topTable3[[i]][(match(results[[i]][,typeID],topTable3[[i]][,typeID])),] = results[[i]]
  }
  return(list(results,topTable3))
}
