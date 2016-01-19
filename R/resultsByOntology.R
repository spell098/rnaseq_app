resultsByOntology = function(selectedOntology,result,ontology){
  resultsOntologies = vector("list",length(selectedOntology))
  names(resultsOntologies) = selectedOntology
  
  for(i in 1:length(selectedOntology)){
    selected <- selectedOntology[i]
    genes <- ontology[[3]][selected]
    resultsOntologies[[i]] <- result[match(genes[[1]],as.character(result$ensembl_gene_id)),]
    resultsOntologies[[i]] <- cbind(resultsOntologies[[i]],rep(selected,nrow(resultsOntologies[[i]])))
    colnames(resultsOntologies[[i]])[ncol(resultsOntologies[[i]])] <- "Ontology"
  }
  resultsOntology = do.call(rbind, resultsOntologies)
  resultsOntology$Ontology = as.character(resultsOntology$Ontology)
  resultsOntology2 = resultsOntology
  for(j in 1:nrow(resultsOntology)){
    duplicates=grep(as.character(resultsOntology[j,]$ensembl_gene_id),as.character(resultsOntology$ensembl_gene_id))
    resultsOntology2[j,]$Ontology = paste(resultsOntology[duplicates,"Ontology"],collapse=",")
  }
  resultsOntology3 = unique(resultsOntology2)
  return(resultsOntology3)
}
