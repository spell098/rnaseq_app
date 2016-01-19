#' Generate boxplots of every groups for a selected gene. Also display stars for significance between groups
#' @author Simon J Pelletier
#' @param selectedGene Gene selected to display boxplots for each groups compared
#' @param names.unique A vector of the names of all groups
#' @param names A vector of the names of all samples
#' @param expr.matrix A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...
#' @param resultsSummary A matrix of data from every gene that is significanlty different
#' between at least two groups. Columns = samples, rows = genes,transcripts,CpG...
#' @param results The results of every comparisons
#' @param color The colors of all boxplot. One color for each group. Could be found in the function, not imported
#' @examples
#' #Illumina HumanHT-12 V4.0 expression beadchip
#' annotation_geo('GPL10558')
#' @keywords geo annotation
#' @seealso
#' \code{\link[GEOquery]{getGEO}} ,
#' \code{\link[GEOquery]{Meta}} ,
#' \code{\link[GEOquery]{Table}} ,
#' \code{\link[annotate]{readGEOAnn}}
#' @note The function should require less input.
#' names and names.unique could be removed by changing the code
#' only the line from the selected gene should be imported, not every object to select it.
#' This would remove expr.matrix,resultsSummary,results,selectedGene
#' colors could also be found in the function. this would leave only a vector from expr.matrix(data) and a vector from resultsSummary(significance)
#' @export
boxplot_element = function(selectedGene,names.unique,names,expr.matrix,resultsSummary,results,color){
  nrows=NULL
  for(i in 1:length(names.unique)){
    nrows[i]=length(grep(names.unique[i],names))
  }
  nrowMax=max(nrows)
  expr=matrix(ncol=length(names.unique),nrow=nrowMax)
  for(i in 1:length(names.unique)){
    expr[1:nrows[i],i] = expr.matrix[selectedGene,grep(names.unique[i],names)]
  }
  colnames(expr) = names.unique
  pvalues = NULL

  height = max(expr[!is.na(expr)])+1
  boxplot(expr,ylim = c(min(expr[!is.na(expr)]),height),col=unique(color))

  coordinates = matrix(ncol=2,nrow=length(results))
  rownames(coordinates) = names(results)
  colnames(coordinates) = c("a","b")
  for(i in 1:nrow(coordinates)){
    coordinates[i,] = match(strsplit(names(results)[i],"-")[[1]],colnames(expr))
  }
  if(nrow(resultsSummary) > 0){
    for(i in 1:(ncol(resultsSummary)-2)){
      pvalues[i] = resultsSummary[match(selectedGene,rownames(resultsSummary)),i+2]
    }
    for(i in 1:nrow(coordinates)){
      if (pvalues[i] < 0.05){
        first = min(coordinates[i,])
        second=max(coordinates[i,])
        ht = max(expr[!is.na(expr)])+0.3+(0.2*(second-first))
        segments(y0=ht,y1=ht,x0=first+0.1,x1=second-0.1,lwd=1)
        if (pvalues[i] < 0.001){
          text(x=mean(c(first,second)), y=ht+0.1, "***", cex=1.2)
        } else if (pvalues[i] < 0.01){
          text(x=mean(c(first,second)), y=ht+0.1, "**", cex=1.2)
        } else text(x=mean(c(first,second)), y=ht+0.1, "*", cex=1.2)
      }
    }
  }
  #barplot(expr,ylim = c(min(expr[!is.na(expr)]),height),col=color)
}
