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
