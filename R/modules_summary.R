modules_summary = function(results,expr.toBind){
  ngenes = nrow(expr.toBind)
  modules.unique = as.numeric(as.character(unique(results$module)))
  modules.all = as.numeric(as.character(expr.toBind[,"module"]))
  modules = as.numeric(as.character(results$module))
  modules_table = data.frame("Module"=1,"Hits"=1,"total"=1,"pvalue"=1,"qvalue"=1,"ratio"=1)
  totalModules = length(unique(expr.toBind[,"module"]))
  for (i in 1:length(modules.unique)){
    n = sum(modules == modules.unique[i])
    N = sum(modules.all == modules.unique[i])
    modules_table[i,c(1,2,3,6)] = c(modules.unique[i],n,N,n/N)
  }
  modules_table = modules_table[order(modules_table[,"Hits"],decreasing = T),]
  rownames(modules_table) = modules_table$Module
  for (c in 1:nrow(modules_table)){         
    # TO CORRECT P-VALUE         
    # # of test = # of annotations tested         
    # p-value of GO is calculated using an hypergeometric distribution         
    hits = modules_table$Hits[c]         
    possibleHits = modules_table$total[c]      
    M = totalModules         
    N = ngenes-possibleHits #total number of genes         
    pvalue = 1-phyper(hits-1,possibleHits,N,nrow(results))                  
    modules_table[,"pvalue"][c] = pvalue
    
  }       
  modules_table = modules_table[order(modules_table[,"pvalue"]),]       
  modules_table[,"qvalue"] = round(p.adjust(modules_table[,"pvalue"],method="BH"),digits=4)
  for (i in 1:length(modules_table[,"qvalue"])){
    if (modules_table[,"qvalue"][i] < 0.001){
      modules_table[,"qvalue"][i] = "< 0.001"
    } 
  } 
  return(modules_table)
}
