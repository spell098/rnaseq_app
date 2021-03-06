% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot_element.R
\name{boxplot_element}
\alias{boxplot_element}
\title{Generate boxplots of every groups for a selected gene. Also display stars for significance between groups}
\usage{
boxplot_element(selectedGene, names, expr.matrix, resultsSummary, results,
  color)
}
\arguments{
\item{selectedGene}{Gene selected to display boxplots for each groups compared}

\item{names}{A vector of the names of all samples}

\item{expr.matrix}{A matrix of standardized data. Columns = samples, rows = genes,transcripts,CpG...}

\item{resultsSummary}{A matrix of data from every gene that is significanlty different
between at least two groups. Columns = samples, rows = genes,transcripts,CpG...}

\item{results}{The results of every comparisons}

\item{color}{The colors of all boxplot. One color for each group. Could be found in the function, not imported}

\item{names.unique}{A vector of the names of all groups}
}
\description{
Generate boxplots of every groups for a selected gene. Also display stars for significance between groups
}
\note{
The function should require less input.
names and names.unique could be removed by changing the code
only the line from the selected gene should be imported, not every object to select it.
This would remove expr.matrix,resultsSummary,results,selectedGene
colors could also be found in the function. this would leave only a vector from expr.matrix(data) and a vector from resultsSummary(significance)
}
\examples{
expr.matrix <- readRDS("data/expr_matrix_LGVD.rds")
results <- readRDS("data/expr_matrix_LGVD.rds")
selectedGene <- "ENSRNOG00000046319"
topTable3 <- readRDS("data/topTable3_LGVD.rds")
resultsSummary <- results_summary(results,topTable3)
color <- c("green","blue","red")
boxplot_element(selectedGene,expr.matrix,resultsSummary,results,color)
}
\author{
Simon J Pelletier
}
\seealso{
\code{\link[GEOquery]{getGEO}} ,
\code{\link[GEOquery]{Meta}} ,
\code{\link[GEOquery]{Table}} ,
\code{\link[annotate]{readGEOAnn}}
}
\keyword{annotation}
\keyword{geo}

