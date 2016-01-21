#' Counts the number of significant transcripts of the same gene
#' @author Simon J Pelletier
#' @import limma
#' @param result A matrix of data. Columns = samples, rows = genes,transcripts,CpG...
#' @param annotations The design of the experiment (object)
#' @param typeID typeID The type of ID used, e.g.
#' @return Linear model contrast ...
#' @keywords limma linear
#' @examples
#' expr.matrix=readRDS("data/expr_matrix_LGVD.rds")
#' results=readRDS("data/results_LGVD.rds")
#' result <- results[[1]]
#' annotations = annotate_ensembl(rownames(expr.matrix))[[1]]
#' transcriptCount <- transcript_count(result,annotations)
#' @export
transcript_count = function(result,annotations,typeID="ensembl_gene_id"){
  count = apply(result,1,function(x){
    transcript_signif = length(grep(x[typeID], result[,typeID]))
    transcript_total = length(grep(x[typeID], annotations[,typeID]))
    ratio_transcript_signif = transcript_signif/transcript_total
    transcript_results = c(transcript_signif,ratio_transcript_signif)
    names(transcript_results) = c("transcript_signif","ratio_transcrpti_signif")
    return(transcript_results)
  })
  colnames(count) = result[,typeID]
  return(t(count))
}
