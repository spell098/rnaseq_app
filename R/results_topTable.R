#' Calculate results
#' @author Simon J Pelletier
#' @import limma
#' @param lm2.contrast MArrayLM object
#' @param expr.toBind Expression data.matrix with supplemental colums: module and ID
#' @param pvalue p-value limit (alpha)
#' @param logFC log2 fold change limit #VERIFY IT IS REALLY LOG2 (default = 1)
#' @param typeID The type of ID used (default="ensembl_gene_id")
#' @param genes_annotation_unique All the annotations for every ID in the 2nd row in expr.toBind
#' (or rownames of expr.matrix)
#' @param adjust Which correction for multiple analysis to use (default = "no").
#' Note: None is different than no somehow
#' @param annotations Annotation (???)
#' @return
#' \describe{
#'  \item{results}{topTable of only the significant results}
#'  \item{topTable3}{Complete topTable}
#' }
#' @keywords limma linear
#' @seealso
#' \code{\link[limma]{topTable}} ,
#' \code{\link[limma]{MArrayLM}} ,
#' \code{\link{transcript_count}} ,
#' \code{\link{expr.toBind}}
#' @examples
#' expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
#' design <- design_contrasts(get_names(colnames(expr.matrix)))
#' results_list <- results_topTable()
#' lm2 <- lm2Contrast(expr.matrix,design)
#' lm2.contrast = lm2[[1]]
#' contrasts=lm2[[2]]
#' contrast.matrix=lm2[[3]]
#' @export
results_topTable = function(lm2.contrast,expr.toBind,pvalue=0.05,logFC=1,typeID="ensembl_gene_id",genes_annotation_unique,annotations,adjust="no"){
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
