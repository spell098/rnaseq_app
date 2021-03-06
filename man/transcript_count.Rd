% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcript_count.R
\name{transcript_count}
\alias{transcript_count}
\title{Counts the number of significant transcripts of the same gene}
\usage{
transcript_count(result, annotations, typeID = "ensembl_gene_id")
}
\arguments{
\item{result}{A matrix of data. Columns = samples, rows = genes,transcripts,CpG...}

\item{annotations}{The design of the experiment (object)}

\item{typeID}{typeID The type of ID used, e.g.}
}
\value{
Linear model contrast ...
}
\description{
Counts the number of significant transcripts of the same gene
}
\examples{
expr.matrix=readRDS("data/expr_matrix_LGVD.rds")
results=readRDS("data/results_LGVD.rds")
result <- results[[1]]
annotations = annotate_ensembl(rownames(expr.matrix))[[1]]
transcriptCount <- transcript_count(result,annotations)
}
\author{
Simon J Pelletier
}
\keyword{limma}
\keyword{linear}

