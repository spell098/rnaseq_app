topTable_unsignificative = function(lm2ReactionContrast,expr.toBind,pvalue,logFC,typeID,specieEnsembl,symbol,selectedComparison){
  
  i = which(colnames(lm2ReactionContrast$coefficients)==selectedComparison)
  toptable = topTable(lm2.contrast , coef=i ,  sort.by = 
                        "none", n=nrow(lm2ReactionContrast$coefficients))
  
  if (!is.null(expr.toBind)) {
    if (colnames(toptable)[1] != "ID"){
      topTable2 = cbind(expr.toBind[,"module"],rownames(toptable),toptable)
    } else {
      topTable2 = cbind(expr.toBind[,"module"],toptable)
    }
    colnames(topTable2)[1:2] = c("module",typeID)  
  } else { 
    if (colnames(toptable)[1] != "ID"){
      topTable2 = cbind(rownames(toptable),toptable)
    } else { 
      topTable2 = toptable
    }
    colnames(topTable2)[1] = typeID   
  }
  
  IDs = biomart_annotation_1(symbol,as.character(topTable2[,typeID]),typeID,specieEnsembl)
  IDs.annotated.unique = IDs[[1]]
  
  resultsReactions = merge(IDs.annotated.unique,topTable2, by=typeID)
  #Pourrait être remplacé par la fonction mean_duplicated
  #mask = !duplicated(results$ensembl_transcript_id)
  #results = results[mask,]
  ratio_transcript_signif = vector("numeric",nrow(resultsReactions))
  transcript_total = vector("numeric",nrow(resultsReactions))
  
  #pas efficace, mais ca marche  # Il y a des replicats de transcripts, pas bon... à fixer
  if(nrow(resultsReactions) >= 1){    
    count_transcript = transcript_count(resultsReactions,IDs,typeID)
    resultsReactions = cbind(resultsReactions,count_transcript)
  }
  return(resultsReactions)
}
