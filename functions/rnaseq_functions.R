affy2ensembl = function(affyID){
  write.csv(affyID,"affyID.csv",quote=FALSE)
  system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$2] == 1 {print $1","$2}\' annotations/affyID.csv annotations/databases/affygene.csv > annotations/affy2ensembl.csv')
  ensemblID = read.csv("affy2ensembl.csv")
  return(ensemblID)
} 

clusteringSummary = function(expr.matrix,pval){
  if(ncol(expr.matrix) < 20 ) ncluster = ncol(expr.matrix) else ncluster = 20
  fuzzy=vector("list",ncluster)
  for(i in 2:ncluster){
    validClusters=NULL
    dominantClusters=NULL
    fuzzy[[i]] = vector("list",4)
    names(fuzzy[[i]]) = c("clusters","clusters_table","correct","correct in subsample")
    clusters = fuzzy_clustering(expr.matrix,i)[[1]]
    clustersUnique = unique(clusters)
    sampleNamesUnique = unique(names(clusters))
    fuzzy[[i]][[1]]=clusters
    fuzzy[[i]][[2]] = as.data.frame(matrix(ncol=length(sampleNamesUnique)*2+4,nrow=length(clustersUnique)+3))
    colnames(fuzzy[[i]][[2]]) = c(sampleNamesUnique,"Dominant samples","% in cluster","% of dominant samples","total in cluster",paste("pvalue",sampleNamesUnique,sep="_"))
    rownames(fuzzy[[i]][[2]]) = c(clustersUnique,"Dominant cluster","% in dominant cluster","total samples")
    
    # j = rows
    # k & l = column
    for (j in 1:length(clustersUnique)){
      for(k in 1:length(sampleNamesUnique)){
        fuzzy[[i]][[2]][j,k]=length(grep(sampleNamesUnique[k],names(clusters[grep(clustersUnique[j],clusters)])))
      }
    }
    for(l in 1:length(sampleNamesUnique)){
      max = max(as.numeric(fuzzy[[i]][[2]][,l][1:length(clustersUnique)]))
      sum = sum(as.numeric(fuzzy[[i]][[2]][,l][1:length(clustersUnique)]))
      fuzzy[[i]][[2]][length(clustersUnique)+1,l] = rownames(fuzzy[[i]][[2]])[as.numeric(fuzzy[[i]][[2]][,l]) == max][1]
      fuzzy[[i]][[2]][length(clustersUnique)+2,l] = max/sum
      fuzzy[[i]][[2]][length(clustersUnique)+3,l] = sum(as.numeric(fuzzy[[i]][[2]][,l][1:length(clustersUnique)]))
    } 
    for(j in 1:length(clustersUnique)){
      
      noCorrect=matrix(ncol=length(sampleNamesUnique),nrow=length(clustersUnique)) #calculate the ratio of the signif group in cluster
      noCorrect2 = noCorrect #calculate the ratio of signif samples of total samples
      colnames(noCorrect) = sampleNamesUnique
      rownames(noCorrect) = clustersUnique
      samplesValues = as.numeric(fuzzy[[i]][[2]][j,][1:length(sampleNamesUnique)])
      totalSamples = as.numeric(fuzzy[[i]][[2]][length(clustersUnique)+3,][1:length(sampleNamesUnique)])
      ratios = samplesValues/totalSamples
      maxRatio=max(ratios)
      sumValues=sum(samplesValues)
      
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+1] = paste(sampleNamesUnique[ratios==maxRatio],collapse=" & ")
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+2] = max(samplesValues)/sumValues
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+3] = max(samplesValues)/sumValues
      fuzzy[[i]][[2]][j,length(sampleNamesUnique)+4] = sumValues
      for(l in 1:length(sampleNamesUnique)){
        #N = nsamples-totalSamples #rest
        m = totalSamples[l]
        k = sumValues
        x=samplesValues[l]
        n=sum(totalSamples) - totalSamples[l]
        pvalue = 1-phyper(x,m,n,k)
        fuzzy[[i]][[2]][j,length(sampleNamesUnique)+4+l] = pvalue
        if(pvalue < pval) dominantClusters[l] = sampleNamesUnique[l]
      }  
    }
    for (z in 1:length(clustersUnique)){
      for (y in 1:length(sampleNamesUnique)){
        pvalue = fuzzy[[i]][[2]][z,length(sampleNamesUnique)+4+y]
        if (pvalue < pval){
          validClusters[z] = sampleNamesUnique[y]
          noCorrect[z,y] = as.numeric(fuzzy[[i]][[2]][z,y])/as.numeric(fuzzy[[i]][[2]][nrow(fuzzy[[i]][[2]]),y])
          noCorrect2[z,y] = as.numeric(fuzzy[[i]][[2]][z,y])/ncol(expr.matrix)
        } else {
          noCorrect[z,y] = 0
          noCorrect2[z,y] = 0
        }
      }
    }
    if(!is.null(validClusters)){
      correctClusters = fuzzy[[i]][[2]][1:length(validClusters),][!is.na(validClusters),]
      totalSamples = NULL
      sampling = rep(0,length(sampleNamesUnique))
      names(sampling) = sampleNamesUnique
      ratioSamples = NULL
      for (z in 1:length(sampleNamesUnique)){
        totalSamples[z] = sum(as.numeric(correctClusters[,z]))
      }
      for (z in 1:length(sampleNamesUnique)){
        for (y in 1:nrow(correctClusters)){
          if(correctClusters[y,(z+length(sampleNamesUnique)+4)] < pval){
            sampling[z] = sampling[z]+as.numeric(correctClusters[y,z])
          }
        }
      }
      fuzzy[[i]][[3]] = sum(noCorrect2)
      ratioSamples = sampling/totalSamples
      names(ratioSamples) = sampleNamesUnique
      fuzzy[[i]][[4]] = ratioSamples
    } else {
      fuzzy[[i]][[4]] = rep(0,length(sampleNamesUnique))
    }
    
  } 
  return(fuzzy)
}

plotCorrectAttribution = function(fuzzy){
  correct = NULL
  for(i in 2:length(fuzzy)){
    correct[i] = fuzzy[[i]][[3]]
  }
  plot(correct,type="o",ylim = c(0,1),axes=FALSE)
  axis(side = 1,at=2*(0:10))
  axis(side = 2,at=0.1*(0:10))
  #lines(correct)
  #text(1:length(fuzz), correct, 1:length(fuzz), cex=0.6, pos=3, col="red")
}
plotCorrectSubsample = function(fuzzy){
  z=length(fuzzy[[4]][[4]])+1
  correct = matrix(nrow=z,ncol=20)
  for(i in 2:length(fuzzy)){
    for(j in 1:length(fuzzy[[i]][[4]])){
      if(!is.null(fuzzy[[i]][[4]][j])){
        correct[j,i] = fuzzy[[i]][[4]][j]
      } else {
        correct[j,i] = 0
      }
    }
  }
  correct[z,] = apply(correct[-(z),],2,function(x){
    #x[1]
    (15*x[1]+35*x[2])/50 # a corriger
  })
  #rownames(correct) = names(fuzzy[[i]][[4]])
  plot(x=NULL,xlim=c(1, 20), ylim=c(0,1))
  for (k in 1:nrow(correct)){
    lines(correct[k,],type="o",col=k)
  }
  layout(rbind(1,2), heights=c(7,1)) 
  #rect(0.5, 0.9-(nrow(correct)*0.02), 2.8, 1)
  legend(x=0,y=1,names(fuzzy[[i]][[4]]),lty=c(1,1),col=1:nrow(correct),cex=0.6,text.width=1)
}

results_summary = function(results,topTable3,adjust,typeID){
  genes = NULL
  if(adjust=="no") sortedBy = "P.Value" else sortedBy = "adj.P.Val"
  for(i in 1:length(results)){
    if(!is.null(results[[i]])){
      genes = unique(c(genes,as.character(results[[i]][,"symbol"])))
    }
  }
  print(genes)
  if(length(genes) > 0){
    if(typeID != "symbol"){ 
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,c("symbol","ensembl_gene_id","entrezgene_id")]
      rownames(resultsSummary) = resultsSummary[,typeID]
    } else {
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,c("symbol","entrezgene_id","ID")]
      rownames(resultsSummary) = resultsSummary[,typeID]
    }
    resultsSummary = resultsSummary[,-2]
    for (i in 1:length(topTable3)){
      resultsSummary[,i+2] = topTable3[[i]][match(genes,topTable3[[i]]$symbol),][,sortedBy]
      colnames(resultsSummary)[i+2] = names(topTable3)[i]
    }
  } else {
    if(typeID != "symbol"){ 
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,c("symbol","ensembl_gene_id","entrezgene_id")]
    } else {
      resultsSummary=topTable3[[1]][match(genes,topTable3[[1]][,"symbol"]),][,typeID]
    }
  }
  return(resultsSummary)
}

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

changeID <- function(expr.matrix,typeID,genes_annotation_unique){
  symbols = as.character(genes_annotation_unique[match(rownames(expr.matrix),genes_annotation_unique[,"ID"]),typeID])
  pos <- !is.na(symbols) & !duplicated(symbols)
  expr.matrix2 <- expr.matrix[pos,]
  rownames(expr.matrix2) = symbols[pos]
  return(expr.matrix2)
}
  
fuzzy_clustering = function(expr.matrix,noCluster){
  mask = !duplicated(rownames(expr.matrix))
  EXPR_MAT_selected= expr.matrix[mask,]
  in_mfuzz = standardise(new("ExpressionSet", exprs = t(as.matrix(EXPR_MAT_selected ))))
  cl_mfuzz <- cmeans(t(EXPR_MAT_selected), noCluster, m=1.25)
  if (length(strsplit(colnames(expr.matrix)[1],"_")[[1]]) > 1){
    names = get_names(colnames(expr.matrix))[[1]]
  } else {
    names = get_names(colnames(expr.matrix))[[1]]
  }
  clusters=cl_mfuzz$cluster
  names(clusters) = names
  return(list(clusters,in_mfuzz,cl_mfuzz))
}

numberNames = function(names){ 
  names=as.character(names)
  names2=NULL
  for (i in unique(names)){
    l=grep(i,names)
    k=0
    for(j in l){
      k=k+1
      names2[j] = paste(names[j],k,sep="_")
    }
  }
  return(names2)
}

volcanoPlot <- function(topTable3,logfc,pval){
  data1 = topTable3
  data2 = data1[!duplicated(data1[,1]) & !duplicated(data1$logFC),]
  signifGenes = subset(data2,P.Value < pval)
  signifGenesAdj = subset(data2,adj.P.Val < pval)
  signifGenesPFC = subset(data2,(abs(logFC) > logfc) & data2$P.Value < pval)
  signifGenesFCAdj = subset(data2,(abs(logFC) > logfc) & data2$adj.P.Val < pval)
  g = ggplot() + 
    geom_point(data=data2, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=1.75) +
    geom_point(data=signifGenes, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=1.75,colour="grey") +
    geom_point(data=signifGenesPFC, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=1.75,colour="blue") +
    geom_point(data=signifGenesAdj, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=3,colour="red") +
    geom_point(data=signifGenesFCAdj, aes(x=logFC, y=-log2(P.Value)), alpha=1, size=3,colour="orange") +
    geom_text(data = signifGenesFCAdj,mapping=aes(x=logFC, y=-log2(P.Value), label=symbol),size=3, vjust=1.5, hjust=0.5) + 
    theme(legend.position = "none") +
    xlab("log2 fold change") + ylab("-log2 p-value")
  g 
  #+ geom_hline(aes(yintercept=-log2(pval))) + geom_vline(aes(xintercept=logFC))
}

MAPlot <- function(topTable3,logfc,pval){
  data1 = topTable3
  data2 = data1[!duplicated(data1[,1]) & !duplicated(data1$logFC),]
  signifGenes = subset(data2,P.Value < pval)
  signifGenesAdj = subset(data2,adj.P.Val < pval)
  signifGenesPFC = subset(data2,(abs(logFC) > logfc) & data2$P.Value < pval)
  signifGenesFCAdj = subset(data2,(abs(logFC) > logfc) & data2$adj.P.Val < pval)
  g = ggplot() + 
    geom_point(data=data2, aes(x=AveExpr, y=logFC), alpha=1, size=1.75) +
    geom_point(data=signifGenes, aes(x=AveExpr, y=logFC), alpha=1, size=1.75,colour="grey") +
    geom_point(data=signifGenesPFC, aes(x=AveExpr, y=logFC), alpha=1, size=1.75,colour="blue") +
    geom_point(data=signifGenesAdj, aes(x=AveExpr, y=logFC), alpha=1, size=3,colour="red") +
    geom_point(data=signifGenesFCAdj, aes(x=AveExpr, y=logFC), alpha=1, size=3,colour="orange") +
    geom_text(data = signifGenesFCAdj,mapping=aes(x=AveExpr, y=logFC, label=symbol),size=3, vjust=1.5, hjust=0.5) + 
    theme(legend.position = "none") +
    ylab("log2 fold change") + xlab("Average Expression")
  g 
  #+ geom_hline(aes(yintercept=-log2(pval))) + geom_vline(aes(xintercept=logFC))
}


PCA = function(expr.matrix,names,colorGroups,namesColor,cex.size){
  rs_rda = rda(t(expr.matrix))
  color=NULL
  for (i in 1:length(unique(names))){
    color[grep(unique(names)[i],names)] = wes()[i]
  }
  plot(rs_rda)
  text(rs_rda, display = "sites", labels = colnames(expr.matrix), choices = c(1, 2), cex= cex.size,col=color,adj=0.5,pos=3)
  legend("topright",legend=unique(namesColor),cex=1,col=unique(colorGroups),lwd=1)
  if(length(color)>length(unique(color)) & sum(color != colorGroups)>0){
    legend("topleft",legend=unique(names),cex=1,col=unique(names),lwd=1)
  }
  points(rs_rda,lwd=1,col=color,bg=colorGroups,pch=21)
  #control = grep("control",colnames(expr.matrix))
  #depression= grep("depression",colnames(expr.matrix))
  #bipolar = grep("bipolar",colnames(expr.matrix))
  #schizophernia = grep("schizophrenia",colnames(expr.matrix))
  #col=c("blue","green","orange","red")
  for(i in 1:length(unique(colnames(expr.matrix)))){
    ordiellipse(rs_rda,groups = colorGroups,show.groups = unique(colorGroups)[i],conf = 0.95,col = unique(colorGroups)[i])
  }
}

RDA_factors <- function(names){
  labs = as.data.frame(do.call("rbind", strsplit(names,"_")))
  colnames(labs) = lapply(strsplit(comparisons,"\\\\"),function(x){
    paste(x,collapse="-")
  })
  mod0 <- rda(t(expr.matrix) ~ 1, labs)
  mod1 <- rda(t(expr.matrix) ~ ., labs)
  #mod <- ordistep(mod0, scope=formula(mod1))
  #plot(mod)
  step.res <- ordiR2step(mod0, mod1, perm.max = 200)
  step.res$anova  # Summary table
  plot(step.res)
}

annotation_geo = function(geoAccession){
  typeID = "symbol" #Rename the symbol column simply "symbol" instead of "rgd_symbol" or "hgnc_symbol"
  annotations = getGEO(geoAccession)
  symbol_column = NULL
  entrez_column = NULL
  iter = 0
  annotation_table <- Table(dataTable(annotations))
  symbol_column <- grep("symbol",colnames(annotation_table),ignore.case=TRUE)
  entrez_column <- grep("entrez",colnames(annotation_table),ignore.case=TRUE)
  iter=iter+1
  if(length(symbol_column) == 0 | length(entrez_column) == 0){
    tmp<-strsplit(Meta(annotations)$relation," ")[[1]]
    geoAccession2<-tmp[grep("GPL",tmp)]
    print("Looking for alternative platform")
    annotation_table<-readGEOAnn(geoAccession2) #Change to make this command only if getGEO return empty table
    symbol_column <- grep("symbol",colnames(annotation_table),ignore.case=TRUE)
    entrez_column <- grep("entrez",colnames(annotation_table),ignore.case=TRUE)
  }
  if(length(symbol_column) == 0 | length(entrez_column) == 0) return(data.frame())
  #symbols <- as.character(annotation_table[,symbol_column])
  colnames(annotation_table)[c(symbol_column,entrez_column)] = c("symbol","entrezgene_id")
  return(annotation_table)
}

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



color_groups = function(names,names.unique){
  colors = vector("character",length(names))
  for (i in 1:length(names.unique)){
    colors[(grep(names.unique[i],names))] = wes()[i] 
  }
  return(colors)
}

WGCNA_modules = function(expr.matrix,ncores,names.unique){
  print("Running WGCNA to find modules")
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads(nThreads=ncores)
  nSets = length(names.unique);
  setLabels = names.unique
  # Form multi-set expression data: columns starting from 9 contain actual expression data.
  multiExpr = vector(mode = "list", length = nSets)
  
  for(i in 1:length(names.unique)){
    multiExpr[[i]] = list(data = as.data.frame(t(expr.matrix[,grep(names.unique[i],colnames(expr.matrix))])))
  }
  exprSize = checkSets(multiExpr)
  powers = c(seq(4,10,by=1), seq(12,20, by=2));
  powerTables = vector(mode = "list", length = nSets);
  # Call the network topology analysis function for each set in turn
  saveRDS(powerTables,file="powerTables.rds")
  
  print("Looking for soft thresholds")
  powerTables = foreach(set=1:nSets) %dopar% {
    #print(paste0("Soft threshold #",set," / ",nSets))
    list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                       verbose = 2)[[2]])
  }
  collectGarbage()
  print("Looking for consensus modules")
  bnet = blockwiseConsensusModules(
    multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
    pamRespectsDendro = FALSE, 
    mergeCutHeight = 0.25, numericLabels = TRUE,
    minKMEtoStay = 0,
    saveTOMs = TRUE, verbose = 5,
    nThreads=ncores)
  save(bnet,file = "bnet.rdata")
  print("Modules informations saved as bnet.rdata")
  return(bnet)
}
colnamesTable = function(x,names){
  colnames(x) <- names
  return(x)
}
geneOntology = function(results,go,typeID,ngenes){
  
  ontologyInfos = vector("list",length(results))
  names(ontologyInfos) = names(results)
  geneGO <- by(go$go_term_name,
               go[,typeID],
               function(x) as.character(x))
  GOgene <- by(go[,typeID],
               go$go_term_name,
               function(x) as.character(x))
  GOgene = GOgene[unique(unlist(geneGO))]
  GOarray = array(dim = length(GOgene))
  
  for (i in 1:length(results)){
    
    if (length(results[[i]]) > 0){
      
      genes=as.character(results[[i]][,typeID])
      significant_geneGO = geneGO[genes]
      significant_geneGO = significant_geneGO[!is.na(names(significant_geneGO))]
      significant_GO = unique(unlist(significant_geneGO))
      significant_GOlist = GOgene[significant_GO]
      GOarray <- lapply(significant_GO,function(x){ 
        names(significant_geneGO[grep(x,significant_geneGO,fixed = TRUE)])
      })
      mask <- significant_GO != ""
      GOarray <- GOarray[mask]
      
      names(GOarray) <- significant_GO[mask]
      GOarray = order_list(GOarray)
      
      a=match(genes,go[,typeID])
      a=a[!is.na(a)]
      GOarray_CellComponent = GOarray[go[a,]$go_domain == "cellular_component"]
      GOarray_MolecularFunction = GOarray[go[a,]$go_domain == "molecular_function"]
      GOarray_BiologicalProcess = GOarray[go[a,]$go_domain == "biological_process"]
      
      
      #Find the most represented GO         
      GO.percent = matrix(ncol = 5,nrow = length(GOarray))
      rownames(GO.percent) = names(GOarray)
      colnames(GO.percent) = c("#Hits","#Possible_genes","p-value","q-value","Hit_ratio")
      
      #Faire un lapply???
      if(length(GOarray) > 0){
        for (c in 1:length(GOarray)){
          # TO CORRECT P-VALUE
          # # of test = # of annotations tested
          # p-value of GO is calculated using an hypergeometric distribution
          hits = length(GOarray[names(GOarray[c])][[1]])
          if(hits > 0){
            possibleHits = length(significant_GOlist[names(GOarray[c])][[1]])
            GO.percent[,"#Hits"][names(GOarray)[c]] = hits
            GO.percent[,"#Possible_genes"][names(GOarray)[c]] = possibleHits
            GO.percent[,"Hit_ratio"][names(GOarray)[c]] = hits/possibleHits
            M = length(significant_GO)
            N = ngenes-possibleHits #total number of genes
            pvalue = 1-phyper(hits-1,possibleHits,N,nrow(results[[i]]))
            if(length(pvalue) == 0) {
              pvalue=1
              print("oops")
            }
            GO.percent[,"p-value"][names(GOarray)[c]] = pvalue  
          } else if(hits <= 0){
            print("PROBLEM WITH #HITS FOUND; SHOULD'NT BE 0")
            GO.percent[,"p-value"][names(GOarray)[c]] = 1  
          }
          
          
        }
        GO.percent = GO.percent[order(GO.percent[,"p-value"]),]
        GO.percent[,"q-value"] = p.adjust(GO.percent[,"p-value"],method="BH")
        
      }
      
      
      #Find the gene with the most functions
      #geneGO.selected
      
      #Create an object with all lists(a list of lists)
      ontologyInfos[[i]] = list(geneGO,GO.percent,GOarray,GOarray_CellComponent,GOarray_MolecularFunction,GOarray_BiologicalProcess)
      names(ontologyInfos[[i]]) = c("geneGO.selected","GOarray.percent","GOarray","GOarray_CellComponent","GOarray_MolecularFunction","GOarray_BiologicalProcess")
    }
  }
  return(ontologyInfos)
}  
export2cytoscape = function(expr.matrix,results,threshold,typeID){
  expr.matrix.selected = expr.matrix[match(as.character(results[,typeID][!is.na(results[,typeID])]),rownames(expr.matrix)),]
  
  expr.adj.selected = adjacency(matrix(as.numeric(t(expr.matrix.selected)),ncol=nrow(expr.matrix.selected),nrow=ncol(expr.matrix.selected)))
  rownames(expr.adj.selected) = rownames(expr.matrix.selected)
  colnames(expr.adj.selected) = rownames(expr.matrix.selected)
  network = exportNetworkToCytoscape(expr.adj.selected, threshold = threshold)
  network$nodeData = merge(network$nodeData,results, by.x = "nodeName", by.y=typeID,all.y=TRUE)
  network$edgeData = cbind(network$edgeData[,1:4],rep("adjacency(WGCNA)",nrow(network$edgeData)))
  colnames(network$edgeData)[ncol(network$edgeData)] = "type"
  
  numChildren = rep(0,nrow(network$nodeData))
  #pourrait peut-etre etre changé pour utiliser with()
  for(n in 1:nrow(network$nodeData)){
    if (!is.na(as.character(network$nodeData[,1][n]))){
      numChildren[n] = length(grep(as.character(network$nodeData[,1][n]),c(as.character(network$edgeData[,1]),as.character(network$edgeData[,2]))))
    } else {
      numChildren[n] = 0
    }
  }
  network$nodeData = cbind(symbol=network$nodeData,numChildren)
  #filter="symbol"
  #write.csv(network$nodeData[,-12],file = paste0("node_rnaseq_",name,".csv"))
  # write All result GO  
  #htmlReport(go[[i]] , file= paste0("ontology_", colnames(contrast.matrix)[i]), ".html")  
  return(network)
}

wes <- function(){
  col1 = wes_palette("Royal1")
  col2 = wes_palette("Royal2")
  col3 = wes_palette("GrandBudapest")
  col4 = wes_palette("Moonrise1")
  col5 = wes_palette("Moonrise2")
  col6 = wes_palette("Moonrise3")
  col7 = wes_palette("Cavalcanti")
  col8 = wes_palette("Chevalier")
  col9 = wes_palette("BottleRocket")
  col10 = wes_palette("Darjeeling")
  col11 = wes_palette("Darjeeling2")
  col = c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11)
  return(col)
}

design_contrasts = function(names){
  #names=strsplit(colnames(expr.matrix),"_")
  
  #names2=unlist(lapply(names,function(x){
  #  paste(x[-length(x)],collapse="_")
  #}))
  #colnames(expr.matrix) = names
  f <- factor(names, levels = unique(names))
  design <- model.matrix(~0+f)
  colnames(design) <- levels(f)
  rownames(design) = names
  return(design)
}

signifRows = function(expr.matrix,resultsSummary){
  expr.matrix.signif = data.matrix(expr.matrix[match(rownames(resultsSummary),rownames(expr.matrix)),])
  colnames(expr.matrix.signif) = colnames(expr.matrix)
  rownames(expr.matrix.signif) = as.character(resultsSummary[,1])
  return(expr.matrix.signif)
}

meanMatrix = function(expr.matrix.signif,names.unique){
  expr.matrix.signif.means = matrix(ncol=length(names.unique),nrow=nrow(expr.matrix.signif))
  
  for (i in 1:length(names.unique)){
    expr.matrix.signif.means[,i] = data.matrix(apply(expr.matrix.signif[,grep(names.unique[i],colnames(expr.matrix.signif))],1,mean))
  }
  colnames(expr.matrix.signif.means) = names.unique
  rownames(expr.matrix.signif.means) = rownames(expr.matrix.signif)
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  return(expr.matrix.signif.means)
}


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

corrExprSetAgilent = function(expr.matrix) {
  expr.matrix = expr.matrix[!is.na(rownames(expr.matrix)),]
  expr.matrix = expr.matrix[grep("^A_",rownames(expr.matrix)),]
  for(i in 1:nrow(expr.matrix)){
    expr.matrix[i,] = apply(matrix(expr.matrix[grep(rownames(expr.matrix)[i],rownames(expr.matrix)),],ncol = ncol(expr.matrix)),2,mean)
  }
  return(expr.matrix)
}

rename_list_data_frame = function(list){
  list <- lapply(list,function(x){ 
    names(x)<-"ensembl_gene_id" 
    x 
  })
  return(list)
}
comparisonsPheno = function(exprset){
  comparisons_table <- pData(exprset)
  comparisons = apply(comparisons_table,2,function(x){
    paste(unique(x),collapse="\\")
  })
  comparisons2 = NULL
  j<-1
  for(i in 1:length(comparisons)){
    if(length(strsplit(comparisons[i],"\\\\")[[1]]) > 1 & length(strsplit(comparisons[i],"\\\\")[[1]]) < ncol(exprs(exprset))){
      comparisons2[j] = comparisons[i]
      names(comparisons2)[j] = names(comparisons)[i]
      j<-j+1
    } else {
      comparisons_table = comparisons_table[,-j]
    }
  }
  comparisons2 = comparisons2[!duplicated(comparisons2)]
  comparisons_table = comparisons_table[!duplicated(lapply(comparisons_table,summary))]
  return(list(comparisons2,comparisons_table))
}
get_names = function(names){
  names1 = rep("",length(names))
  for (i in 1:length(names)){
    nameSplit = strsplit(names[i],"_")
    for (j in 1:(length(nameSplit[[1]]))){
      if(j==1) names1[i] = nameSplit[[1]][j]
      else if (is.na(as.numeric(gsub("([0-9]+).*$", "\\1", nameSplit[[1]][j])))) names1[i] = paste0(names1[i],"_",nameSplit[[1]][j])
    }
  }
  names.unique = unique(names1)
  return(list(names1,names.unique))
}

all_entrez = function(specieEnsembl){
  ensembl=useMart("ensembl")
  ensembl = useDataset(specieEnsembl,mart=ensembl)
  all_entrez = getBM(attributes=c("ensembl_gene_id","entrezgene"),
                     mart = ensembl)
  all_entrez = all_entrez[!is.na(all_entrez[,2]),][,2]
  return(all_entrez)
}

biomart_annotation_cpg = function(cpg){
  if(length(cpg)==2){
    cpg=as.list(cbind(cpg,cpg[,2]))
  } else {
    cpg=as.list(cpg)
  }
  chr=""
  for(i in 1:length(cpg[,1])){
    chr[i]=strsplit(as.character(cpg[,1][i]),"chr")[[1]][2]
  }
  position=list(as.numeric(chr),cpg[,2],cpg[,3])
  ensembl=useMart("ensembl")
  ensembl = useDataset(specieEnsembl,mart=ensembl)
  ID.annotated = getBM(attributes=c("ensembl_gene_id","chromosome_name",
                                    "start_position","end_position"),
                       filters = c("chromosome_name","start","end"),
                       #filters = "efg_agilent_wholegenome_4x44k_v1",
                       values = position,
                       mart = ensembl)
  
}

heatmap_expr = function(results,expr.matrix.annotated){
  col_clu=redgreen(100)
  col_clu = col_clu[100:1]
  
  signifGenes = unique(c(results[[1]]$symbol,results[[2]]$symbol,results[[3]]$symbol))
  expr.matrix.annotated.signif = expr.matrix.annotated[expr.matrix.annotated$symbol  %in% signifGenes,]
  expr.matrix.signif = data.matrix(expr.matrix.annotated.signif[,2:12])
  expr.matrix.signif.means = cbind(apply(expr.matrix.signif[,1:4],1,mean),apply(expr.matrix.signif[,5:8],1,mean),apply(expr.matrix.signif[,9:11],1,mean))
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  heatmap.2(expr.matrix.signif.means, col=col_clu, scale="row",
            main = "Heatmap result", key=TRUE, symkey=FALSE, 
            density.info="none", trace="none", cexRow=0.5,cex=0.5, 
            dendrogram="column", Colv = TRUE, Rowv = FALSE)
  
}


#Pourrait être beaucoup plus efficace...
mean_duplicated = function(){
  for (i in 1:length(expr.matrix.annotated$ensembl_transcript_id)){
    x = grep(expr.matrix.annotated$ensembl_transcript_id[i],expr.matrix.annotated$ensembl_transcript_id)
    selected = as.matrix(expr.matrix.annotated[x,][,1:ncol(expr$A)+1])
    selected = matrix(as.numeric(selected),ncol = ncol(expr$A))
    mean = apply(selected,2,mean)
    mean.rbind = rbind(mean,mean,mean)
    expr.matrix.annotated[x,1:ncol(expr$A)+1] = mean
  }
  return(expr$A)
}

rename_list_data_frame = function(list){
  list <- lapply(list,function(x){ 
    names(x)<-"ensembl_gene_id" 
    x 
  })
  return(list)
}

merge_list_data_frames_with_data_frame = function(list,results){
  list <- lapply(list,function(x){
    merge(x,results,by = "ensembl_gene_id")
  })
  return(list)
}

order_list = function(list){
  list = list[order(sapply(list,length),decreasing=T)]
  return(list)
}

expr.matrix.annotated.means = function(results,expr.matrix.annotated){
  col_clu=redgreen(100)
  col_clu = col_clu[100:1]
  
  signifGenes = unique(c(results[[1]]$symbol,results[[2]]$symbol,results[[3]]$symbol))
  signifGenes = signifGenes[signifGenes != ""]
  expr.matrix.annotated.signif = expr.matrix.annotated[match(signifGenes,expr.matrix.annotated$symbol),]
  rownames(expr.matrix.annotated.signif) = expr.matrix.annotated.signif$symbol
  expr.matrix.signif = data.matrix(expr.matrix.annotated.signif[,3:13])
  expr.matrix.signif.means = cbind(apply(expr.matrix.signif[,1:4],1,mean),apply(expr.matrix.signif[,5:8],1,mean),apply(expr.matrix.signif[,9:11],1,mean))
  colnames(expr.matrix.signif.means) = names(results)
  #heatmap.2(data.matrix(expr.matrix.annotated.signif[,2:12]))
  return(expr.matrix.signif.means)
}

transcript_count = function(result,annotations,typeID){
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

#ordiellipse(rs_rda,groups = names,show.groups = "DEP_DBS",conf = 0.95,col = "blue")
#ordiellipse(expr.rda,groups = groups,show.groups = "DG_DIETE",conf = 0.95,col = "green")
#ordiellipse(expr.rda,groups = groups,show.groups = "DG_INSULIN",conf = 0.95,col = "red")

#labs      =  do.call("rbind", strsplit(colnames(expr.matrix),"_"))
#maternal_care   = labs[,1]
#maternal_care   = unlist(lapply(maternal_care,as.character)) 
#brain_part= labs[,2]
#brain_part [brain_part == "NA"] <- "A"
#brain_part= unlist(lapply(brain_part,as.character)) 
#rda.vegan = rda(t(expr.matrix) ~ maternal_care + brain_part )
#rda.vegan
#rs_vegan_matrix = get.result.vegan(rda.vegan)
#plot(rda.vegan,choices = c(1, 2), scaling=1, display = c("sp","cn"), cex = 0.2, col = "white", main= "PC1_PC2_scl1")

#Dendogramme
#I should do some bootstrap thing to have most probable groups
#matrix_distance = dist(t(expr.matrix),method="euclidian")
#rs.hcl          = hclust(matrix_distance, method="average")
#library(cclust)
#plot(rs.hcl, sub="",labels = NULL, main = "Dendrogramme",xlab = "labels", ylab = "Hauteur", cex= 0.5, las = 2)
#t.expr <- t(expr.matrix)
#cl = cclust( t.expr, 3, method="kmeans")
#plot(matrix(apply( t.expr,1,mean),nrow=8,ncol=1),type="n",xlim=c(0,10))
#text(matrix(apply( t.expr)),labels=rownames(t.expr),cex=0.5,col=cl$cluster)

#cluster analysis

lm2Contrast <- function(expr.matrix,design){
  lm <- lmFit(expr.matrix, design)
  colnames(design) = as.character(colnames(design))
  contrasts = ""
  n = 1
  for (i in 1:(length(colnames(design))-1)){
    for (j in (i+1):length(colnames(design))){
      contrasts[n] = paste0(colnames(design)[i],"-",colnames(design)[j])
      n = n + 1
    }
  }
  contrast.matrix <-makeContrasts(contrasts=contrasts, levels = design)
  lm2 = contrasts.fit(lm, contrast.matrix)
  lm2.contrast <- eBayes(lm2)
  
  return(list(lm2.contrast,contrasts,contrast.matrix))
}

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

exprToBind = function(bnet,expr.matrix,typeID){
  expr.toBind = cbind(bnet$colors,rownames(expr.matrix),expr.matrix)
  colnames(expr.toBind)[1:2] = c("module",typeID)
  return(expr.toBind)
}

expr.matrix.annotation = function(ID.annotated.unique,expr.toBind){
  colnames(expr.toBind)[1:2] = c("module","agilentID")
  expr.matrix.annotated = merge(expr.toBind,ID.annotated.unique,by.x = "agilentID", by.y = "efg_agilent_wholegenome_4x44k_v1")
  return(expr.matrix.annotated)
}


cpg.annotation = function(cpg,annotation){
  list_cpg=vector("list",nrow(cpg))
  for (i in 1:nrow(cpg)){
    for (j in 1:nrow(annotation)){
      if ( (cpg[i,2]<annotation[j,3] && cpg[i,3]>annotation[j,2]) || (cpg[i,2]>annotation[j,3] && cpg[i,3]<annotation[j,2]) && as.character(cpg[i,1]) == as.character(annotation[j,1])){
        test=1
        if(is.null(list_cpg[[i]])){
          n=1
        }
        list_cpg[[i]][n] = as.character(annotation[j,4])
        n=n+1
      }
    }
  }
  return(list_cpg)
}


updateCytoscapeNetwork = function(network,nodes,edges){
  network$nodeData = merge(network$nodeData,nodes,by = "symbol",all.x=TRUE)
  network$edgeData = merge(network$edgeData,edges,by.y = c("src","dest"),by.x = c("fromNode","toNode"),all.x=TRUE)
}

updateCytoscapeAll = function(network,nodes,edges){
  network$nodeData = merge(network$nodeData,nodes,by = "symbol",all=TRUE)
  network$edgeData = merge(network$edgeData,edges,by.y = c("src","dest"),by.x = c("fromNode","toNode"),all=TRUE)
  return(network)
}

createReactionCytoscape = function(network,nodes,edges,topTable3){
  network$nodeData = merge(topTable3,nodes,by = "symbol",all.y=TRUE)
  network$edgeData = merge(network$edgeData,edges,by.y = c("src","dest"),by.x = c("fromNode","toNode"),all.y=TRUE)
  return(network)
}

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

get_comparisons = function(names){
  comparisons=NULL
  conditions = matrix(nrow=length(names),ncol = length(strsplit(names[1],"_")[[1]]))
  
  for (i in 1:length(names)){
    conditions[i,] = strsplit(names[i],"_")[[1]]
  }
  for (j in 1:ncol(conditions)){
    tmp=unique(conditions[,j])
    comparisons[j] = paste(tmp,collapse="\\")
    #iter=0
    #for (k in 1:(length(tmp)-1)){
    #  rest=tmp[(k+1):length(tmp)]
    #  for(r in rest){
    #    iter=iter+1
    #    comparisons[iter] = paste(tmp[k],r,sep="\\")
    #  }
    #}
  }
  
  return(comparisons)
}

possibleComparisons = function(names){
  possibleComparisons = NULL
  n=1
  for (i in 1:(length(unique(names))-1)){
    for (j in (i+1):length(unique(names))){
      possibleComparisons[n] = paste0(unique(names)[i],"-",unique(names)[j])
      n=n+1
    }
  }
  return(possibleComparisons)
}

conditionsChoice = function(comparisonSelection,names){ 
  if (length(comparisonSelection) > 0){
    names2=NULL
    for (m in length(comparisonSelection):1){
      conditionsSelected = strsplit(comparisonSelection[m],'\\\\')
      for (j in conditionsSelected[[1]]){
        range=grep(j,names)
        names2[range] = j
      }
    }
  }
  return(names2)
}
namesSelectedComparisons = function(namesTable){
  namesSelected <- apply(namesTable,1,function(x){
    y<-lapply(x,function(z){
      paste(strsplit(z," ")[[1]],collapse="_")
    })
    paste(y,collapse="_")
  })
  names(namesSelected) = rownames(namesTable)
  return(namesSelected)
}



network.graphite = function(pSymbol,expr.matrix,toptableNodes,typeID,threshold,selectedReaction,toptable3,selectedComparison,expr.toBind,pvalue){
  nodes = cbind(symbol=nodes(pSymbol),Reaction=rep(selectedReaction,length(nodes(pSymbol))))
  colnames(nodes) = c("symbol","reaction")
  
  edgeColor = NULL
  arrow=NULL
  tail=NULL
  toto=NULL
  networkSelectedReaction = export2cytoscape(expr.matrix,toptableNodes,threshold,typeID)
  networkSelectedReaction$edgeData[,"fromNode"] = as.character(toptableNodes[,"symbol"])[match(networkSelectedReaction$edgeData[,"fromNode"],toptableNodes[,"ensembl_gene_id"])]
  networkSelectedReaction$edgeData[,"toNode"] = as.character(toptableNodes[,"symbol"])[match(networkSelectedReaction$edgeData[,"toNode"],toptableNodes[,"ensembl_gene_id"])]
  
  edges = merge(networkSelectedReaction$edgeData,edges(pSymbol),by.x=c("fromNode","toNode","type","direction"),by.y=c("src","dest","type","direction"),all.x=TRUE,all.y=TRUE)
  
  gSymbol <- pathwayGraph(pSymbol,edge.type=NULL)
  g2 <- addEdge(networkSelectedReaction$edgeData[,"fromNode"], networkSelectedReaction$edgeData[,"toNode"], gSymbol, networkSelectedReaction$edgeData[,"weight"])
  
  g3 <- layoutGraph(g2)
  
  toptableNodes = merge(toptable3[[selectedComparison]],nodes,by="symbol",all.y=TRUE)
  toptableNodes$logFC[is.na(toptableNodes$logFC)] = 0
  fc=toptableNodes$logFC
  toptableNodes$adj.P.Val[is.na(toptableNodes$adj.P.Val)] = 1
  #modules=rep(1,nrow(toptableNodes))
  if (!is.null(expr.toBind)) {
    toptableNodes$module[is.na(toptableNodes$module)] = 0
    rainbow=rainbow(max(as.numeric(as.character(toptableNodes$module)))+1)
    modules=apply(toptableNodes,1,function(x){
      rainbow[(as.numeric(as.character(x["module"])))+1]
    })
  } else modules = rep("black",nrow(toptableNodes))
  textCol = unlist(modules)
  colorIntensity=apply(toptableNodes,1,function(x){
    if (as.numeric(x["logFC"]) > 0) {col=rgb(red=1,green=0,blue=0,as.numeric(x["logFC"])/(max(toptable3[[selectedComparison]][,"logFC"])+0.000001))
    }else col=rgb(red=0,green=0,blue=1,(abs(as.numeric(x["logFC"])/min(toptable3[[selectedComparison]][,"logFC"])+0.000001)))
  })
  signif=apply(toptableNodes,1,function(x){
    if (as.numeric(x["adj.P.Val"]) < pvalue){ signif="triangle" 
    }else signif="circle"
  })
  
  types.edges = as.character(unique(edges[,"type"]))
  #types.colors = wes[1:length(types.edges)]
  types.colors = c("blue","red","black","cyan","green","dodgerblue","forestgreen","deeppink")
  names(types.colors) = types.edges
  numChildren = networkSelectedReaction$nodeData$numChildren
  names(numChildren) = networkSelectedReaction$nodeData[,4]
  
  
  
  for (i in 1:nrow(edges)){
    a=rownames(edges)[!is.na(match(edges[,2],edges[,1][i]))]
    b=rownames(edges)[!is.na(match(edges[,1],edges[,2][i]))]
    c=match(a,b)
    d=c[!is.na(c)]
    d=a[!is.na(a[c])]
    
    if (length(d) > 0 && as.character(edges[,"type"][i]) != "binding"){
      if(length(grep("ACTIVATION",as.character(edges[i,"type"]))) > 0 && length(grep("ACTIVATION",as.character(edges[d,"type"]))) > 0){ 
        arrow[i] = "normal"
        tail[i]="normal"
        edgeColor[i] = "green"
      } else if (length(grep("INHIBITION",as.character(edges[i,"type"]))) > 0 && length(grep("INHIBITION",as.character(edges[d,"type"]))) > 0){ 
        arrow[i] = "tee"
        tail[i]="tee"
        edgeColor[i] = "green"
      } else if (length(grep("INHIBITION",as.character(edges[i,"type"]))) > 0 && length(grep("ACTIVATION",as.character(edges[d,"type"]))) > 0){ 
        arrow[i] = "tee"
        tail[i]="normal"
        edgeColor[i] = "green"
      } else if (length(grep("ACTIVATION",as.character(edges[i,"type"]))) > 0 && length(grep("INHIBITION",as.character(edges[d,"type"]))) > 0){ 
        arrow[i] = "tee"
        tail[i]="normal"
        edgeColor[i] = "green"
      } else {
        tail[i]="normal"
        arrow[i] = "normal"
        edgeColor[i] = "pink"
      }
        
    }else if(length(grep("ACTIVATION", as.character(edges$type[i]))) == 1){
      arrow[i] = "normal"
      tail[i] = "none"   
      edgeColor[i] = "black"
    } else if (length(grep("INHIBITION", as.character(edges$type[i]))) == 1){
      arrow[i] = "none"
      tail[i] = "tee"
      edgeColor[i] = "cyan"
    } else if (edges$type[i] == "binding"){
      arrow[i] = "none"
      tail[i] = "none"
      edgeColor[i] = "red"
    } else if (edges$type[i] == "adjacency(WGCNA)"){
      arrow[i] = "none"
      tail[i] = "none"
      edgeColor[i] = "blue"
    } else {
      arrow[i] = "none"
      tail[i] = "none"
      edgeColor[i] = "yellow"
    }
  }
  
  
  edgeWeight=apply(edges,1,function(x){
    if (!is.na(x["weight"])) as.numeric(x["weight"])*3 else 1
  })
  edgeNames = paste0(edges[,"fromNode"],"~",edges[,"toNode"])
  names(edgeWeight) = edgeNames
  names(edgeColor) = edgeNames
  names(tail) = edgeNames
  names(arrow) = edgeNames
  
  
  as.numeric(toptable3[[selectedComparison]][,"logFC"])
  names(colorIntensity) = as.character(toptableNodes[,1])
  names(signif) = as.character(toptableNodes[,1])
  names(textCol) = as.character(toptableNodes[,1])
  
  nodeRenderInfo(g3) <- list(
    fill     = colorIntensity,
    col      = "black",
    lty      = "solid",
    lwd      = (numChildren)+1,
    shape    = signif,
    fontsize = 15,
    textCol=textCol,
    cex=0.9)
  
  
  edgeRenderInfo(g3) <-  list(
    col = edgeColor,
    lwd = edgeWeight,
    arrowhead=arrow,
    arrowtail=tail)
  graph.par(list(graph=list(
    main     = selectedComparison,
    sub      = selectedReaction, 
    cex.main = 1.8,
    cex.sub  = 1.4, 
    col.sub  = "gray")
    
  )
  )
  return(g3)
}
selectTopReaction = function(selectedComparison,topReactions){
  x = strsplit(selectedComparison,"-")[[1]]
  alt1=paste(x[1],x[2],sep="-")
  alt2=paste(x[2],x[1],sep="-")
  if (!is.null(topReactions[[alt1]])){
    selectedTopReactions = topReactions[[alt1]]
  } else { 
    selectedTopReactions = topReactions[[alt2]]
  }
  return(selectedTopReactions)
}

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

colorNumericValues <- function(values){
  values2 = values + abs(min(values))
  values3 = values2/max(values2)
  f <- colorRamp(c("green","white", "red"))
  colors <- rgb(f(values3)/255)
  return(colors)
}

addColorsMatrix <- function(resultsContrast,selectedComparison){
  resultsContrastRamp = data.matrix(resultsContrast[,selectedComparison])
  colnames(resultsContrastRamp) = selectedComparison
  rownames(resultsContrastRamp) = rownames(resultsContrast)
  resultsContrastRamp[,selectedComparison] = colorNumericValues(resultsContrast[,comparison])
  return(resultsContrastRamp)
}
addColorsMatrixResults <- function(results,selectedComparison){
  resultsContrastRamp = data.matrix(results[[selectedComparison]]$logFC)
  colnames(resultsContrastRamp) = selectedComparison
  rownames(resultsContrastRamp) = as.character(results[[selectedComparison]]$symbol)
  resultsContrastRamp[,selectedComparison] = colorNumericValues(resultsContrastRamp)
  return(resultsContrastRamp)
  
}


annotate_ensembl = function(IDs){
  print("Homemade annotation")
  write.csv(IDs,"annotations/gene_names.csv",quote=FALSE)
  write.table(IDs,"annotations/gene_names.ssv",quote=FALSE,sep=";")
  #write.table(rownames(expr.matrix),"annotations/gene_names.txt",sep="\t",quote=FALSE)
  if (length(grep("ENSRNOG",IDs[1])) > 0){
    typeID = "ensembl_gene_id"
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1","$3","$4}\' annotations/gene_names.csv annotations/databases/rat_symbols.csv > annotations/genes_annotations.csv')
    system('awk -F ";" \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1","$3","$4}\' annotations/gene_names.ssv annotations/databases/allGO.ssv > annotations/go.csv')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) = c(typeID,"symbol","entrezgene_id")
    go=read.csv("annotations/go.csv",header=FALSE)
    colnames(go) = c(typeID,"go_term_name","go_domain")
  }else if (length(grep("ENSRNOT",IDs[1])) > 0){
    typeID = "ensembl_transcript_id"
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$2] == 1  && ( $3 || $4 ) {print $1","$2","$3","$4}\' annotations/gene_names.csv annotations/databases/rat_symbols.csv > annotations/genes_annotations.csv')
    system('awk -F ";"\'NR==FNR{c[$2]=1;next} c[$2] == 1  && ( $3 || $4 ) {print $1","$2","$3","$4}\' annotations/gene_names.ssv annotations/databases/allGO.ssv > annotations/go.csv')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) = c("ensembl_gene_id",typeID,"symbol","entrezgene_id")
    go=read.csv("annotations/go.csv",header=FALSE)
    colnames(go) = c("ensembl_gene_id",typeID,"go_term_name","go_domain")
  }else if (length(grep("A_",IDs[1])) > 0){
    typeID = "efg_agilent_wholegenome_4x44k_v1"
    #REMPLACER LES ROWNAMES PAR NOMS TRANSCRIPTS
    system("annotations/annotation_agilent_go.sh",input=rownames(expr.matrix))
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
  } else if (length(grep("ENSG",IDs[1])) > 0){
    typeID = "ensembl_gene_id"
    system('awk -F "," \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1","$3","$4}\' annotations/gene_names.csv annotations/databases/hgnc.csv > annotations/genes_annotations.csv')
    system('awk -F ";" \'NR==FNR{c[$2]=1;next} c[$1] == 1  && ( $3 || $4 ) {print $1","$3","$4}\' annotations/gene_names.ssv annotations/databases/go_human.ssv > annotations/go.csv')
    annotations=read.csv("annotations/genes_annotations.csv",header=FALSE)
    colnames(annotations) = c(typeID,"symbol","entrezgene_id")
    go=read.csv("annotations/go.csv",header=FALSE)
    colnames(go) = c(typeID,"go_term_name","go_domain")
  }
  genes_annotation_unique = annotations[!duplicated(as.character(annotations[,typeID])),]
  return(list(annotations,go,genes_annotation_unique,typeID))
}
