calculateBvalues = function(methy,names){
  bvalues = matrix(ncol = ((ncol(methy)-4)/3),nrow=nrow(methy))
  for(i in 1:((ncol(methy)-4)/3)){
    x=((i*3)+3)
    bvalues[,i] <- cbind(as.numeric(methy[,x])/(as.numeric(methy[,(x-1)])+100))
  }
  colnames(bvalues) = names
  rownames(bvalues) = paste0(methy$chr,".",methy$start)
  bvalues[is.na(bvalues)] <- 0
  return(bvalues)
}

findCpgByGene = function(cpg,symbol){
  cpgByGene=NULL
  for(x in as.character(cpg[,symbol])){
    tmp    = grep(x,cpg[,symbol])
    number = length(tmp)
    cpgs   = paste(as.character(cpg[tmp,]$cpg_pos),collapse=",")
    cpgByGene = rbind(cpgByGene,c(x,number,cpgs))
  }
  colnames(cpgByGene) = c(symbol,"#CpGs","cpg position")
  cpgByGene = unique(as.data.frame(cpgByGene))
  return(cpgByGene)
}

methylDifferentialAnalysis = function(meth){
  myDiffList = vector("list",round(nrow(meth)/20000))
  myDiffList=foreach(i=1:round(nrow(meth)/20000)) %dopar% {
    print(paste0("Iteration ",i, " of ",round(nrow(meth)/20000)))
    if (i==round(nrow(meth)/20000)){
      print("last iteration")
      diffMethyl=calculateDiffMeth(meth[(20000*(i-1)+1):nrow(meth),])
    } else {
      diffMethyl=calculateDiffMeth(meth[(20000*(i-1)+1):(20000*i),])
    }

  }
  myDiff=data.frame()
  for (x in myDiffList){
    myDiff = rbind(myDiff,x)
  }
  return(myDiff)
}

methylDifferentialAnalysisAll = function(meth,pc,names){
  myDiff2 = list()
  for(x in pc){
    print("Calculating differential analysis...")
    split = strsplit(x,"-")[[1]]
    attributes(meth)$treatment[grep(split[1],names)] = 1
    attributes(meth)$treatment[grep(split[2],names)] = 0
    myDiff2[[x]] = methylDifferentialAnalysis(meth)
  }
  return(myDiff2)
}

cpg_annotation = function(cpg,transcripts_coord,transcripts_symbols,typeID){
  ### VERIFY genome version (rn5 vs rn6)
  
  cpg_diff_annotated = list()
  list_cpg_outside_gene = list()
  for(x in names(cpg)){
    #CpGs that are in more than 1 transcript are found more than once and are meant to be kept
    list_cpg=cpg_gene_awk(cpg[[x]])
    list_cpg = cbind(list_cpg,paste0(list_cpg$chr_cpg,".",cpg_pos = c(list_cpg$pos_cpg)))
    colnames(list_cpg)[length(list_cpg)] = "cpg_pos"
    
    
    cpg_diff=merge(list_cpg,cpg[[x]],by="cpg_pos")[,c(-9:-12)]
    cpg_diff$ensembl_transcript_id = unlist(lapply(strsplit(as.character(cpg_diff$ensembl_transcript_id),"_"),function(x){
      x[1]
    }))
    
    list_cpg_outside_gene[[x]] = cpg[[x]][is.na(!match(cpg[[x]]$cpg_pos,cpg_diff$cpg_pos)),]
    print(1)
    cpg_diff_coord=merge(cpg_diff,transcripts_coord,by="ensembl_transcript_id")
    print(2)
    cpg_diff_annotated[[x]] = merge(cpg_diff_coord,transcripts_symbols,by="ensembl_transcript_id")
  }
  return(list(cpg_diff_annotated,list_cpg_outside_gene))
}


bind_cpg_pos = function(myDiff,methDiffNeg,methDiffPos,qvalue){
  cpg2 = list()
  for (x in names(myDiff)){
    mask = (as.numeric(myDiff[[x]][,"meth.diff"]) <= methDiffNeg | as.numeric(myDiff[[x]][,"meth.diff"]) >= methDiffPos) & myDiff[[x]]$qvalue <= qvalue
    cpg = myDiff[[x]][mask,]
    cpg=cbind(cpg,apply(cpg,1,function(x) {
      tmp=matrix(paste0(as.character(x[1]),".",as.character(x[2]),collapse=""),ncol=1)
      gsub("[[:space:]]", "", tmp)
    }))
    colnames(cpg)[ncol(cpg)]=("cpg_pos")
    cpg$cpg_pos=as.character(cpg$cpg_pos)
    cpg2[[x]] = cpg
  }
  return(cpg2)
}

queryCpg = function(allCpgs,selectedPos,selectedComparison,transcripts_coord,transcripts_symbols,typeID){
  selectedCpg = NULL
  selectedCpg[["query"]] = allCpgs[[selectedComparison]][match(selectedPos,allCpgs[[selectedComparison]]$cpg_pos),]
  annotated_cpg = cpg_annotation(selectedCpg,transcripts_coord,transcripts_symbols,typeID)
  return(annotated_cpg)
}

queryGene = function(allCpgs,selectedGeneCpg,selectedComparison,transcripts_symbols,typeID){
  selectedGenes=transcripts_symbols[match(selectedGeneCpg,transcripts_symbols[,typeID]),]
  write.table(allCpgs[[selectedComparison]],"annotations/allCpgs.txt",quote=FALSE)
  write.table(selectedGenes,"annotations/selected_genes.txt",quote=FALSE)
  
  system('awk \'BEGIN{i=0}NR==FNR{i++;a[i]=$2;b[i]=$3;c[i]=$5;next }; {for (j = 1; j <= i; j++) if(a[j] == $2 && b[j] > $3 && b[j] < $4) print a[j]"."b[j]"\t"a[j]"\t"b[j]"\t"$3"\t"$4"\t"$5}\' annotations/allCpgs.txt annotations/selected_genes.txt > annotations/selected_genes_cpg.txt')
  
  annotations_cpg=read.table("annotations/selected_genes_cpg.txt",sep="\t",header=FALSE)
  colnames(annotations_cpg) = c("cpg_pos","chr","pos","transcript_start","transcript_end","ensembl_transcript_id")
  
  return(annotations_cpg)
}

cpg_gene_awk = function(cpg2){
  write.csv(cpg2,"annotations/signif_cpgs.csv",quote=FALSE)       
  system('awk -F "," \'BEGIN{i=0}NR==FNR{i++;a[i]=$2;b[i]=$3;next }; {for (j = 1; j <= i; j++) if(a[j] == $1 && b[j] > $2 && b[j] < $3) print a[j]"\t"b[j]"\t"$4"\t"$1"\t"$2"\t"$3"\tpromoter(3000)"}\'  annotations/signif_cpgs.csv annotations/databases/rn5_annotation_promoters_upstream_3000.txt > annotations/cpg_gene.txt')
  system('awk -F "," \'BEGIN{i=0}NR==FNR{i++;a[i]=$2;b[i]=$3;next }; {for (j = 1; j <= i; j++) if(a[j] == $1 && b[j] > $2 && b[j] < $3) print a[j]"\t"b[j]"\t"$4"\t"$1"\t"$2"\t"$3"\tcoding_exon"}\'  annotations/signif_cpgs.csv annotations/databases/rn5_annotation_coding_exons.txt >> annotations/cpg_gene.txt')
  system('awk -F "," \'BEGIN{i=0}NR==FNR{i++;a[i]=$2;b[i]=$3;next }; {for (j = 1; j <= i; j++) if(a[j] == $1 && b[j] > $2 && b[j] < $3) print a[j]"\t"b[j]"\t"$4"\t"$1"\t"$2"\t"$3"\t5-UTR"}\'  annotations/signif_cpgs.csv annotations/databases/rn5_annotation_5utr.txt >> annotations/cpg_gene.txt')
  system('awk -F "," \'BEGIN{i=0}NR==FNR{i++;a[i]=$2;b[i]=$3;next }; {for (j = 1; j <= i; j++) if(a[j] == $1 && b[j] > $2 && b[j] < $3) print a[j]"\t"b[j]"\t"$4"\t"$1"\t"$2"\t"$3"\t3-UTR"}\'  annotations/signif_cpgs.csv  annotations/databases/rn5_annotation_3utr.txt >> annotations/cpg_gene.txt')
  system('awk -F "," \'BEGIN{i=0}NR==FNR{i++;a[i]=$2;b[i]=$3;next }; {for (j = 1; j <= i; j++) if(a[j] == $1 && b[j] > $2 && b[j] < $3) print a[j]"\t"b[j]"\t"$4"\t"$1"\t"$2"\t"$3"\tintron"}\'  annotations/signif_cpgs.csv annotations/databases/rn5_annotation_introns.txt >> annotations/cpg_gene.txt')
  annotations_cpg=read.table("annotations/cpg_gene.txt",sep="\t",header=FALSE)
  colnames(annotations_cpg) = c("chr_cpg","pos_cpg","ensembl_transcript_id","part_chr","part_start","part_end","part")
  return(annotations_cpg)
}

cpg.gene = function(cpg){
  #Trouver TOUS les cpgs pour un gène (pour avoir %methylation)
  #ajouter: 
  #1: promotor : selon https://www.genomatix.de/online_help/help_gems/FAQ_answers.html, 500bp downstream, 100bp upstream (even if exon?)
  #2: exon/intron
  #3: intergene
  

  annotation=read.table("rn5/rn5_annotation.txt")
  annotation_exons=read.table("rn5/rn5_annotation_coding_exons.txt")
  annotation_introns=read.table("rn5/rn5_annotation_introns.txt")
  annotation_3utr=read.table("rn5/rn5_annotation_3utr.txt")
  annotation_5utr=read.table("rn5/rn5_annotation_5utr.txt")
  annotation_promoters_3000=read.table("rn5/rn5_annotation_promoters_upstream_3000.txt")
  
  list_cpg=vector("list",nrow(cpg))
  foreach(i=1:nrow(cpg)) %dopar% {
    print(paste0("Gene #",i," of ",nrow(cpg)))
    list_cpg[[i]]=as.data.frame(matrix(ncol=9))
    n=1
    for (j in 1:nrow(annotation_promoters_3000)){
      if ((cpg[i,2]<annotation_promoters_3000[j,3] && cpg[i,2]>annotation_promoters_3000[j,2]) && as.character(cpg[i,1]) == as.character(annotation_promoters_3000[j,1])){
        list_cpg[[i]][n,1] = strsplit(as.character(annotation_promoters_3000[j,4]),"_")[[1]][1]
        list_cpg[[i]][n,2] = as.character(cpg[i,1])
        list_cpg[[i]][n,3] = paste0(as.character(cpg[i,1]),".",cpg[i,2])
        list_cpg[[i]][n,4] = annotation_promoters_3000[j,2]        
        list_cpg[[i]][n,5] = annotation_promoters_3000[j,3]
        list_cpg[[i]][n,6] = as.character(annotation_promoters_3000[j,6])   
        list_cpg[[i]][n,7] = "promoters(-3000)"
        list_cpg[[i]][n,8] = j 
        n=n+1
      }
    }
    for (k in 1:nrow(annotation)){
      if ((cpg[i,2]<annotation[k,3] && cpg[i,2]>annotation[k,2]) && as.character(cpg[i,1]) == as.character(annotation[k,1])){
        for (j in 1:nrow(annotation_3utr)){
          if ((cpg[i,2]<annotation_3utr[j,3] && cpg[i,2]>annotation_3utr[j,2]) && as.character(cpg[i,1]) == as.character(annotation_3utr[j,1])){
            list_cpg[[i]][n,1] = strsplit(as.character(annotation_3utr[j,4]),"_")[[1]][1]
            list_cpg[[i]][n,2] = as.character(cpg[i,1])
            list_cpg[[i]][n,3] = paste0(as.character(cpg[i,1]),".",cpg[i,2])
            list_cpg[[i]][n,4] = annotation_3utr[j,2]        
            list_cpg[[i]][n,5] = annotation_3utr[j,3]
            list_cpg[[i]][n,6] = as.character(annotation_3utr[j,6])   
            list_cpg[[i]][n,7] = "3-UTR"
            list_cpg[[i]][n,8] = j
            n=n+1
          }
        }
        for (j in 1:nrow(annotation_5utr)){
          if ((cpg[i,2]<annotation_5utr[j,3] && cpg[i,2]>annotation_5utr[j,2]) && as.character(cpg[i,1]) == as.character(annotation_5utr[j,1])){
            list_cpg[[i]][n,1] = strsplit(as.character(annotation_5utr[j,4]),"_")[[1]][1]
            list_cpg[[i]][n,2] = as.character(cpg[i,1])
            list_cpg[[i]][n,3] = paste0(as.character(cpg[i,1]),".",cpg[i,2])
            list_cpg[[i]][n,4] = annotation_5utr[j,2]        
            list_cpg[[i]][n,5] = annotation_5utr[j,3]
            list_cpg[[i]][n,6] = as.character(annotation_5utr[j,6])   
            list_cpg[[i]][n,7] = "5-UTR"
            list_cpg[[i]][n,8] = j
            n=n+1
          }
        }
        for (j in 1:nrow(annotation_exons)){
          if ((cpg[i,2]<annotation_exons[j,3] && cpg[i,2]>annotation_exons[j,2]) && as.character(cpg[i,1]) == as.character(annotation_exons[j,1])){
            list_cpg[[i]][n,1] = strsplit(as.character(annotation_exons[j,4]),"_")[[1]][1]
            list_cpg[[i]][n,2] = as.character(cpg[i,1])
            list_cpg[[i]][n,3] = paste0(as.character(cpg[i,1]),".",cpg[i,2])
            list_cpg[[i]][n,4] = annotation_exons[j,2]        
            list_cpg[[i]][n,5] = annotation_exons[j,3]
            list_cpg[[i]][n,6] = as.character(annotation_exons[j,6])   
            list_cpg[[i]][n,7] = "exon"
            list_cpg[[i]][n,8] = j
            n=n+1
          } 
        }
        if (is.null(list_cpg[[i]][n,7]) || is.na(list_cpg[[i]][n,7])){
          list_cpg[[i]][n,1] = as.character(annotation[k,4])
          list_cpg[[i]][n,2] = as.character(cpg[i,1])
          list_cpg[[i]][n,3] = paste0(as.character(cpg[i,1]),".",cpg[i,2])
          list_cpg[[i]][n,4] = annotation[k,2]        
          list_cpg[[i]][n,5] = annotation[k,3]
          list_cpg[[i]][n,6] = as.character(annotation[k,6])   
          list_cpg[[i]][n,7] = "intron"
          list_cpg[[i]][n,8] = k
          n=n+1
        }   
      }
    }
    print(list_cpg[[i]])
  }
  return(list_cpg)
}




cpg.gene.p = function(cpg,annotation){
  list_cpg=vector("list",nrow(cpg))
  foreach(i=1:nrow(cpg)) %dopar% {
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

#Make a list of all genes that have different CpGs, with:
#The number of different cpg
#Total cpg
#Vector containing all cpgs position
#Vector containing all cpgs positions of different cpgs

genes_cpgs = function(list_cpg,cpg){
  diffMethGenes=unique(unlist(list_cpg))
  positions_diff=vector("list",length(diffMethGenes))
  positions_all=vector("list",length(diffMethGenes))
  for(i in 1:length(diffMethGenes)){
    #positions_all[[diffMethGenes[i]]]=myobj[grep(diffMethGenes[i],list_cpg),2]
    
    positions_diff[[diffMethGenes[i]]]=cpg[grep(diffMethGenes[i],list_cpg),2]
  }
  return(positions_diff)
}

plot_genes_cpg = function(positions_diff,genes,meth){
  
  #Faut chercher pour le bon chromosome aussi
  #
  if (length(genes) == 1) {
    data=getData(meth)[match(positions_diff[[genes]],getData(meth)[,2]),]
  } else {
    data=getData(meth)[match(genes,getData(meth)[,2]),]
  }
  data_rel=cbind(data[,6]/data[,5],data[,9]/data[,8],data[,12]/data[,11],data[,15]/data[,14])
  rownames(data_rel)=1:nrow(data)
  colnames(data_rel) = c("dorsal","dorsal","ventral","ventral")
  
  data_group=data.frame()
  for(i in 1:nrow(data_rel)){
    data_group=rbind(data_group,cbind(data_rel[i,],paste0(data[i,1],".",data[i,2]),colnames(data_rel)))
  }
  colnames(data_group)=c("value","CpG","group")
  rownames(data_group)=1:nrow(data_group)
  data_group[,1] = as.numeric(levels(data_group[,1]))[data_group[,1]]
  data_group[,2]=as.character(data_group[,2])
  
  #data_group[,2] = as.numeric(levels(data_group[,2]))[data_group[,2]]
  data_group[,3]=as.character(data_group[,3])
  
  boxplot(value~CpG*group,data=data_group,col=rainbow(nrow(data_rel)),las=2,cex.axis=0.5)
  
  library(ggplot2)
  library(ggplot2)
  
  data.summary=data.frame()
  n=1
  for(i in 1:length(unique(data_group[,"CpG"]))){
    for(j in 1:length(unique(data_group[,"group"]))){
      CpG=unique(data_group[,"CpG"])
      group=unique(data_group[,"group"])
      data.summary[n,"CpG"]=CpG[i]
      data.summary[n,"group"]=group[j]
      popo=data_group[which(data_group[,"CpG"]==CpG[i]),]  
      data.summary[n,"mean"]=mean(popo[which(popo[,"group"]==group[j]),1])
      data.summary[n,"sd"]=sd(popo[which(popo[,"group"]==group[j]),1])
      data.summary[n,"se"]=data.summary[n,"sd"]/sqrt(2)
      #data.summary[n,"ci"]=confint(popo[which(popo[,"group"]==group[j]),1])
      n=n+1
    }
  }
  data.summary$group <- factor(data.summary$group)
  data.summary$CpG <- factor(data.summary$CpG)
  
  
  ggplot(data.summary, aes(x=CpG, y=mean, fill=group)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}

data_genes_all_positions = function(annotation,genes,methy){
  ok = rep(0,nrow(positions))
  annotation_gene=cbind(annotation[match(genes,annotation[,4]),],genes)
  positions=annotation_gene[,c(1:3,13)]  
  list_positions=vector("list",nrow(positions))
  names(list_positions)=as.character(positions[,4])
  #Find all CpGs
  for (i in 1:nrow(positions)){
    n=1
    print(paste0("Gene #",i," (",positions[i,4],") of ",nrow(positions)))
    for (j in 1:nrow(methy)){
      if (((methy[j,2]<positions[i,3] && methy[j,3]>positions[i,2]) || 
             (methy[j,3]>positions[i,2] && methy[j,2]<positions[i,3])) &&
            as.character(methy[j,1]) == as.character(positions[i,1])){
        list_positions[[as.character(positions[i,4])]][n] = paste0(methy[j,1],".",methy[j,2])
        n=n+1
      } 
    }
  }
}

all.cpg.gene = function(methy,annotation,gene){
  genes_positions=annotation[match(gene,annotation[,4]),]
  list_cpg=vector("list",nrow(genes_positions))
  for (i in 1:nrow(genes_positions)){
    print(paste0("Gene #",i," of ",nrow(genes_positions)))
    for (j in 1:nrow(annotation)){
      if ( (genes_positions[i,2]<methy[j,3] && genes_positions[i,3]>methy[j,2]) || (genes_positions[i,2]>methy[j,3] && genes_positions[i,3]<methy[j,2]) && as.character(genes_positions[i,1]) == as.character(methy[j,1])){
        if(is.null(list_cpg[[i]])){
          n=1
        }
        list_cpg[[i]][n] = as.numeric(levels(methy[j,3]))[methy[j,3]]
        n=n+1
      }
    }
  }  
}

bvalues.signif = function(bvalues,cpg_annotated,names,names.unique,typeID,symbol){
  bvalues.signif = NULL
  for (i in 1:length(cpg_annotated)){
    bvalues.signif[[i]] = bvalues[match(cpg_annotated[[i]]$cpg_pos,rownames(bvalues)),]
    rownames(bvalues.signif[[i]]) = cpg_annotated[[i]][,symbol]
    bvalues.signif[[i]] = bvalues.signif[[i]][rownames(bvalues.signif[[i]]) != "",]
    #rownames(expr.matrix.annotated.signif) = cpg_annotated[[i]]$rgd_symbol[match(rownames(expr.matrix.signif),cpg_annotated[[i]]$ensembl_gene_id)]
  } 
  
  bvalues.signif2 = unique(do.call(rbind, bvalues.signif))
  bvalues.signif.means = matrix(ncol=length(names.unique),nrow=nrow(bvalues.signif2))
  
  for (i in 1:length(names.unique)){
    bvalues.signif.means[,i] = data.matrix(apply(bvalues.signif2[,grep(names.unique[i],colnames(expr.matrix))],1,mean))
  }
  colnames(bvalues.signif.means) = names.unique
  rownames(bvalues.signif.means) = rownames(bvalues.signif2)
  return(list(bvalues.signif2,bvalues.signif.means))
}

