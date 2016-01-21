#  title: "rnaseqs_analysis"
#author: Simon J. Pelletier
#Status  Public on Jul 01, 2011
#Title:
#  High-frequency stimulation of the ventrolateral thalamus regulates gene expression in hippocampus, motor cortex and caudate-putamen
#Organism  Rattus norvegicus
#Experiment type  Expression profiling by array
#Summary  :
#  Transcriptional profiling of rat hippocampus comparing rats subjected to chronic (1 hour daily during 14 days) vetrolateral thalamus stimulation (DBS group) with sham operated rats (with electrode placement but without stimulation) and with naive rats.

#Overall design  Three-condition experiment, DBS, sham and naive. Every sample consisted of pooled hippocampi from three rats. Slices between Bregma levels -2 and -4 were used for dissections of the hippocampi. Common reference design, comparing every condition versus a brain pooled tissue consisting of pooled hippocampal, thalamic and cortical brain tissue of the naive (n=4), sham-operated (n=5) and stimulated (n=4) rats. Biological replicates: 4 DBS, 4 sham and 3 naive. One replicate per array.
#Deep brain stimulation
#Notes:
#Il faudrait voir si pour trouver les networks si on devrait faire des sous-groupes selon ce qu'on veut étudier: exemple, high_dorsal - low_dorsal : on regarderait les gènes qui varient ensemble dans seulement ce sous-groupe

library(rnaseqApp)
library(foreach)
library(parallel)
library(doParallel)
library(doMC)
library(foreign)
library(limma)
library(gdata)
library(biomaRt)
library(AnnotationHub)
library(GenomicRanges)
library(vegan)
library(gplots)
library(GOstats)
library(pathview)
library(WGCNA)
#library(rgb)
library(wesanderson)
library(cclust)
library(RCurl)
library(parallel)
library(plyr)
library(graphite)
library(SPIA)
library(RCytoscape)
library(org.Rn.eg.db)
library(Rgraphviz)
library(network)
ncores = detectCores()
cl = makeForkCluster(nnodes = ncores)
registerDoParallel(cl)
registerDoMC()

dataFileName = "expr.matrix.VvsD.rdata"
specieEnsembl = "rnorvegicus_gene_ensembl"
symbol = "rgd_symbol"
organism = "rnorvegicus"

wes=wes()
expr.matrix=readRDS("data/expr_matrix_LGVD.rds")
#expr.matrix=readRDS("expr_matrix.rds")
expr.matrix = expr.matrix[!duplicated(rownames(expr.matrix)),]
annotation1 = annotate_ensembl(rownames(expr.matrix))

annotations=annotation1[[1]]
go=annotation1[[2]]
genes_annotation_unique = annotation1[[3]]
typeID=annotation1[[4]]
names = get_names(expr.matrix)
comparisons = get_comparisons(names)
comparisonSelection = comparisons[1]
#expr.matrix = conditionsChoice(comparisonSelection,expr.matrix,names)
#possibleComparisons(names)

######   Normalization ################
#print("Normalization of the data")
design = design_contrasts(expr.matrix)
#voom2 = voom(expr.matrix,design, plot=F, normalize="quantile")
#voom3 = voom(expr.matrix,design, plot=T)
#expr.matrix = removeBatchEffect(as.matrix(voom2$E), batch=c("1","2","1","2","1","2","1","2"), design=voom2$design)
#saveRDS(expr.matrix,"data/rnaseq_exprMatrix.rds")
#expr.matrix=voom2
#######################################
#bnet=WGCNA_modules(expr.matrix)
#save(bnet,file="rnaseq_bnet.rdata")
#saveRDS(bnet,"rnaseq_bnet.rds")

bnet<-readRDS("data/bnet_LGVD.rds")
expr.toBind = cbind(bnet$colors,rownames(expr.matrix),expr.matrix)
colnames(expr.toBind)[1:2] = c("module",typeID)

print("Running linear model")
lm2=lm2Contrast(expr.matrix,design)
#saveRDS(lm2,file="rnaseq_lm2.rdata")
lm2.contrast = lm2[[1]]
contrasts=lm2[[2]]
contrast.matrix=lm2[[3]]

print("Finding significative results")
logFC=c(-1.3,1.3)
pvalue=0.05
results_list = results_topTable(lm2.contrast,expr.toBind,pvalue,logFC,typeID,genes_annotation_unique,annotations,"no")
results = results_list[[1]]
topTable3 = results_list[[2]]
saveRDS(results,file="data/results_LGVD.rds")
saveRDS(topTable3,file="data/topTable3_LGVD.rds")

#rownames(lm2.contrast$coef) = ID.annotated.unique[,1][match(rownames(lm2.contrast$coef),ID.annotated.unique[,4])]

print("Finding gene ontologies")

ontology1 = geneOntology(results,go,typeID,nrow(expr.matrix))
saveRDS(ontology,file="data/rnaseq_ontology.rds")


selectedOntology = c("membrane","protein","cell","nucleus")
selectedComparison = "HIGH_DORSAL-HIGH_VENTRAL"
resultsOntology = resultsByOntology(selectedOntology,results[[selectedComparison]],ontology1)
saveRDS(resultsOntology,file="data/rnaseq_resultsOntology.rdata")

modules_table = modules_summary(results[[selectedComparison]],expr.toBind)


selectedModule = 28
maskModule=as.numeric(results[[selectedComparison]]$module) == selectedModule
resultsModule = results[[selectedComparison]][maskModule,]
ontologyModule = geneOntology(resultsModule,go)


threshold = 0.7
filter = "rgd_symbol"
name=names(results)[1]
network = export2cytoscape(expr.matrix,results[[selectedComparison]],threshold,typeID)

write.csv(network$nodeData,file = paste0("node_rnaseq_",selectedComparison,".csv"))
write.csv(network$edgeData,file = paste0("edge_rnaseq_",selectedComparison,".csv"))

networkSelectedOntology = export2cytoscape(expr.matrix,resultsOntology,threshold,filter,typeID)
write.csv(networkSelectedOntology$nodeData,file = paste0("node_rnaseq_",paste0(selectedComparison,"_",selectedOntology),".csv"))
write.csv(networkSelectedOntology$edgeData,file = paste0("edge_rnaseq_",paste0(selectedComparison,"_",selectedOntology),".csv"))



# kegg, biocarta, reactome
#reactions = reactions(results,ratReactome,specieEnsembl)

#saveRDS(reactions,paste0("hippocampusVD","_reactions.rds"))
topReactions = readRDS(paste0("data/hippocampusVDHL","_reactions.rds"))

ratReactome <- pathways("rnorvegicus", "reactome")
print("preparing for SPIA...")
prepareSPIA(reactome, "reactome",print.names=TRUE)

selectedReaction = "Degradation of the extracellular matrix"
p <- ratReactome[[selectedReaction]]
pSymbol <- convertIdentifiers(p, "SYMBOL")
toptableNodes = cbind(topTable3[[selectedComparison]][match(nodes(pSymbol),topTable3[[selectedComparison]]$rgd_symbol),],Reaction = rep(selectedReaction,length(nodes(pSymbol))))
toptableNodes = toptableNodes[!is.na(toptableNodes[,1]),]
network=network.graphite(pSymbol,expr.matrix,toptableNodes,typeID,threshold,selectedReaction,topTable3,selectedComparison)
renderGraph(network)
legend(-300,900,legend = c("KO","WT"),cex=0.5, lwd=c(2.5,2.5),col =c("pink","green"),bty="n",x.intersp=0.3)

networkReaction = createReactionCytoscape(network,nodes,edges,topTable3[[selectedComparison]])
network
#networkReaction = updateCytoscape(network,nodes,edges)
#IL FAUDRA QUE LES NODES DES REACTIONS AIENT LES DATA, LISTE DE RESULTATS UPDATÉE
updateResults
#networkReaction = updateCytoscapeAll(network,nodes,edges)
updateNodesData = function(resultsReactions)

