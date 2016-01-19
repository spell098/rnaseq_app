reactions = function(results,entrezgene_ids,specieEnsembl){
  ALL <- as.character(entrezgene_ids[!is.na(entrezgene_ids)])
  print(head(ALL))
  res=vector("list",length(results))
  names(res) = names(results)
  print("Looking for top reactions...")
  for (i in 1:length(results)){
    print(paste("running results list",i,"of",length(results),sep=" "))
    if (nrow(results[[i]]) > 0){
      DE = as.numeric(results[[i]][!is.na(results[[i]][,"entrezgene_id"]),]$logFC)
      names(DE) <- as.character(as.vector(results[[i]][!is.na(results[[i]]$entrezgene),]$entrezgene))
      print(DE)
      x = match(ALL,names(DE))
      DE = DE[x[!is.na(x)]]
      print(DE)
      if(!is.null(DE)){
        res[[i]] <- runSPIA(de=DE, all=ALL, "reactome")
        print(paste0("DONE #",i))
      } else {
        print("NULL")
      }
    } else {
      res[[i]] = NULL
      print(paste0("DONE #",i,"\nNothing Found..."))
    }
  }
  print("All top reactions finished")
  return(res)
}