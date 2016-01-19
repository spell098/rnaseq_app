# I D'ONT USE THIS FUNCTION ANYMORE
all_entrez = function(specieEnsembl){
  ensembl=useMart("ensembl")
  ensembl = useDataset(specieEnsembl,mart=ensembl)
  all_entrez = getBM(attributes=c("ensembl_gene_id","entrezgene"),
                     mart = ensembl)
  all_entrez = all_entrez[!is.na(all_entrez[,2]),][,2]
  return(all_entrez)
}
