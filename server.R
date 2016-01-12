
# All large files are on ratking and retrieved using the system(wget ...,intern=TRUE) function 


#setwd("~/workbench/codes/r/phd/rnaseq_analysis2/rnaseq_analysis")
#library(annotate)
#library(XMLRPC)
library(GEOquery)
library(foreach)
#library(DT)
#library(methylKit)
#library(parallel)
library(doParallel)
#library(doMC) 
#library(foreign)
library(limma)
library(gdata)
#library(biomaRt)
library(Mfuzz)
#library(GenomicRanges)
library(vegan)
library(gplots)
library(GOstats)
library(pathview)
source("rnaseq_functions.R")
source("methylAnalysis_functions.R")
library(WGCNA)
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
library(shiny)
library(a4Base)
library(IRanges)
library(graphite)
require(ggplot2)

#library(rsconnect)

#WHEN THE DESIGN IS CHANGED (GROUPS CHANGED) THE SAMPLES NEED TO BE RENORMALIZED, 
#THE PIPELINE HAS TO START FROM RAW DATA TO ALLOW THAT


shinyServer(
  function(input, output, session) {
    #ncores = detectCores()
    #cl = makeForkCluster(nnodes = ncores)
    #registerDoParallel(cl)
    #registerDoMC()
    
    ################################################################################
    ###                              RNASEQ                                      ###
    ################################################################################
    
    col_clu=redgreen(100)
    col_clu = col_clu[100:1]    
    
    observeEvent(input$submit1, {
      
      topReactions = reactive({NULL})
      modulesTable = reactive({NULL})
      expr.toBind = reactive({NULL})
      summaryClust = reactive({NULL})
      
      specieEnsembl=reactive({input$specieEnsembl})
      if(input$fileContent == "Expression Matrix"){
        expr.matrix = reactive({readRDS(input$File1[,"datapath"])})
        names1 = reactive({get_names(colnames(expr.matrix()))[[1]]})
        comparisons = reactive({get_comparisons(names1())})
        names.unique = reactive({get_names(colnames(expr.matrix()))[[2]]})
        annotation1 = reactive({annotate_ensembl(rownames(expr.matrix()))})
        annotations=reactive(annotation1()[[1]])
        go=reactive({annotation1()[[2]]})
        genes_annotation_unique = reactive({annotation1()[[3]]})
        typeID=reactive({annotation1()[[4]]})
        clusters1=reactive({NULL})
      } else if(input$fileContent == "Expression Set"){
        if(length(input$geo_id) > 0  & input$submit1 == 1){
          gset <- reactive({getGEO(input$geo_id, GSEMatrix =TRUE)}) #GSE61276 GSE12654
          exprset <- reactive({gset()[[1]]})
        } else if (length(input$geo_id) == 0){
          exprset <- reactive({readRDS(input$File1[,"datapath"])})
        }
        expr.matrix <- reactive({exprs(exprset())})
        names1 <- reactive({get_names(sampleNames(exprset()))[[1]]})
        names.unique <- reactive({get_names(sampleNames(exprset()))[[2]]})
        comparisons <- reactive({comparisonsPheno(exprset())[[1]]})
        comparisonsTable <- reactive({comparisonsPheno(exprset())[[2]]})
        annotation1 <- reactive({annotation_geo(annotation(exprset()))})
        goGPL <- reactive({annotation1()[,c(colnames(annotation1())[grep("ontology",colnames(annotation1()),ignore.case=TRUE)],"symbol")]})
        observe({
          write.table(goGPL(),"goTable.txt",sep=";",quote=FALSE,col.names=FALSE,row.names=FALSE)
        })
        annotations <- reactive({annotation1()[,c("ID","symbol","entrezgene_id")]})
        system("annotateGPL.sh")
        go<-reactive({readRDS("goTable3.txt")})
        typeID <- reactive({"symbol"})
        genes_annotation_unique <- reactive({annotations()[!duplicated(as.character(annotations()[,typeID()])),]})
      }
      
      if (input$groupBy == "fuzzy"){
        clusters1 = reactive({fuzzy_clustering(expr.matrix(),input$noCluster)[[1]]})
        color = reactive({unlist(lapply(clusters1(),function(x){wes()[x]}))})
        namesColor = reactive({paste0("a",clusters1())})
        comparisons = reactive({get_comparisons(namesColor())})
        clusterNames = reactive({numberNames(namesColor())})
        selectedVariables = reactive({comparisons()[1]})
        names2 = eventReactive(input$submit1,{conditionsChoice(selectedVariables(),namesColor())})
      } else {
        if(input$fileContent == "Expression Matrix"){
          if(input$submit1 <= 1){
            names2=reactive({names1()})
            namesColor = reactive({names1()})
            clusterNames = reactive({names2()})
          } else {
            names2 = eventReactive(input$submit1,{conditionsChoice(input$selectedVariables,names1())})
            namesColor = eventReactive(input$submit1,{names1()})
            clusterNames = reactive({names2()})
          }
        } else if(input$fileContent == "Expression Set"){
          if(input$submit1 > 1){
            namesTable <- reactive({comparisonsTable()[,input$selectedVariables]})
            selectedVariables = reactive({input$selectedVariables})
          } else {
            namesTable <- reactive({comparisonsTable()})
            selectedVariables = reactive({colnames(namesTable())})
          }
          names2 <- reactive({namesSelectedComparisons(namesTable())})
          #names2 <- namesSelectedComparisons(namesTable)
          #names2 = reactive({conditionsChoice(selectedVariables(),namesSelected())})
          namesColor = eventReactive(input$submit1,{names2()})
          clusterNames = reactive({names2()})
        }
      }
      names.unique2 = reactive({unique(names2())})
      design = eventReactive(input$submit1,{design_contrasts(names2())})
      labs = reactive({do.call("rbind", strsplit(names2(),"_"))})
      
      
      if (!is.null(input$File2)){
        bnet = reactive({readRDS(input$File2[,"datapath"])})
        expr.toBind = reactive({exprToBind(bnet(),expr.matrix(),typeID())})
      } else {
        expr.toBind = reactive({cbind(expr.matrix(),module=rep("white",nrow(expr.matrix())))})
      }
      
      
      observeEvent(input$submitMethylModule,{
        bnet = reactive({WGCNA_modules(expr.matrix())})
        expr.toBind = reactive({exprToBind(bnet(),expr.matrix(),typeID())})
      })
      
      lm2=reactive({lm2Contrast(expr.matrix(),design())})
      lm2.contrast = reactive({lm2()[[1]]})
      contrasts=reactive({lm2()[[2]]})
      contrast.matrix=reactive({lm2()[[3]]})
      
      #if(input$submit1 == 0){
      #  if(input$fileContent == "Expression Matrix"){
      #    possibleComp = reactive({possibleComparisons(names1())})
      #  } else if(input$fileContent == "Expression Set"){
      #    possibleComp = reactive({possibleComparisons(namesSelectedComparisons(namesTable()))})
      #  }
      #} else {
      possibleComp = eventReactive(input$submit1,{possibleComparisons(names2())})
      #}
      
      summaryClust = reactive({clusteringSummary(expr.matrix(),input$pval_clust)})
      
      selectedComparison = reactive({input$selectedComparison})
      pvalue = eventReactive(input$submit2,{input$pvalue})
      logFC = eventReactive(input$submit2,{c(input$logFCneg,input$logFCpos)})
      results_list = reactive({results_topTable(lm2.contrast(),expr.toBind(),pvalue(),logFC(),typeID(),genes_annotation_unique(),annotations(),input$adjust)})
      results = reactive({results_list()[[1]]})
      topTable3 = reactive({results_list()[[2]]})
      resultsSummary = reactive({results_summary(results(),topTable3(),input$adjust,typeID())})
      RowSideColors=reactive({expr.toBind()[match(rownames(resultsSummary()),rownames(expr.toBind())) ,"module"]})
      if(input$groupBy == "groups"){
        color = reactive({color_groups(rownames(lm2.contrast()$design),colnames(lm2.contrast()$design))})
      }  
      
      expr.matrix.signif = reactive({signifRows(expr.matrix(),resultsSummary())})
      expr.matrix.signif.means = reactive({meanMatrix(expr.matrix.signif(),names.unique2())})
      
      resultsContrast = reactive({signifRows(lm2.contrast()$coefficients,resultsSummary())})
      
      #if (is.null(input$topGO_rows_selected) == TRUE){
      #  selectedOntology = reactive({strsplit(input$selectedOntology,"\n")[[1]]})
      #} else {
      selectedOntology = reactive({c(input$topGO_rows_selected,strsplit(input$selectedOntology,"\n")[[1]])})
      #}
      
      
      ontology1 = reactive({geneOntology(results(),go(),typeID(),nrow(genes_annotation_unique()))})
      topGO = reactive({cbind(rownames(ontology1()[[selectedComparison()]][[2]]),ontology1()[[selectedComparison()]][[2]])})
      resultsOntology = reactive({resultsByOntology(selectedOntology(),results()[[selectedComparison()]],ontology1()[[selectedComparison()]])})
      
      networkSelectedComparison = reactive({export2cytoscape(expr.matrix(),results()[[selectedComparison()]],input$threshold,typeID())})
      networkSelectedOntology = reactive({export2cytoscape(expr.matrix(),resultsOntology(),input$threshold,typeID())})
      
      if (!is.null(input$File2)){
        modulesTable = reactive({modules_summary(results()[[selectedComparison()]],expr.toBind())})
        maskModule=reactive({as.numeric(as.character(results()[[selectedComparison()]]$module)) == input$modulesTable_rows_selected})
        resultsModule=reactive({NULL})
        #selectedModuleOntology = reactive({input$topGO_rows_selected})
        resultsModule = reactive({list(resultsModule=results()[[selectedComparison()]][maskModule(),])})
        ontologyModule = reactive({geneOntology(resultsModule(),go(),typeID(),length(rownames(expr.matrix())))})
        topModuleGO = reactive({cbind(rownames(ontologyModule()[[1]][[2]]),ontologyModule()[[1]][[2]])})
        resultsModuleOntology = reactive({resultsByOntology(input$topModuleGO_rows_selected,resultsModule()[[1]],ontologyModule()[[1]])})

      }
      
      if (!is.null(input$File3)){
        topReactions = reactive({readRDS(input$File3[,"datapath"])})
        selectedTopReactions = reactive({selectTopReaction(selectedComparison(),topReactions())})
      }
      ratReactome <- reactive({pathways(input$specie, "reactome")})
      observeEvent(input$findTopReactions,{
        print("preparing for SPIA...")
        prepareSPIA(ratReactome(), "reactome",print.names=TRUE)
        print("Looking for top reactions...")
        topReactions = reactive({reactions(results(),unique(annotations()$entrezgene_id),specieEnsembl())})
        saveRDS(topReactions(),paste0("data/reactions_",paste(input$selectedVariables,collapse = "_"), ".rds"))
      })    
      
      selectedReaction <- reactive({c(strsplit(input$selectedReactions,"\n")[[1]],strsplit(selectedTopReactions()[input$topReactions_rows_selected,1],"\n")[[1]])})
      p <- reactive({ratReactome()[[selectedReaction()]]})
      pSymbol <- reactive({convertIdentifiers(p(), "SYMBOL")})
      toptableNodes = reactive({
        cbind(
          topTable3()[[selectedComparison()]][match(nodes(pSymbol()),topTable3()[[selectedComparison()]]$symbol),],
          Reaction = rep(selectedReaction()[[1]],length(nodes(pSymbol())))
        )
      })
      toptableNodes1 = reactive({toptableNodes()[!is.na(toptableNodes()[,1]),]})
      
      #toptableNodes = reactive({cbind(toptableNodes(),Reaction = rep(selectedReaction(),nrow(toptableNodes())))})
      
      network1=reactive({network.graphite(pSymbol(),expr.matrix(),toptableNodes1(),typeID(),input$threshold,selectedReaction()[[1]],topTable3(),selectedComparison(),NULL,pvalue())})
      
      ################################################################################
      ###                              OUTPUTS                                     ###
      ################################################################################
      output$text1 <- renderText({ 
        print((clusters1()))
        #print(selectedTopReactions()[[selectedComparison()]][input$topReactions_rows_selected,])
      })
      ################################################################################
      ###                              PLOTS                                       ###
      ################################################################################
      output$boxplot_element <- renderPlot({
        if (is.null(input$resultsSummary_rows_selected) == TRUE){
          x=input$selectedGene
        } else {
          x=input$resultsSummary_rows_selected
        }
        boxplot_element(x,names.unique2(),names2(),expr.matrix(),resultsSummary(),results(),color())
      })
      output$boxplot_element_raw <- renderPlot({
       x=input$selectedGenes
        boxplot_element(x,names.unique2(),expr.matrix(),resultsSummary(),topTable3(),color())
      })
      output$pathway <- renderPlot({
        renderGraph(network1())
      },height = 800, width = 1200)
      
      
      output$venn <- renderPlot({vennDiagram(resultsContrast(),cex=input$venn.cex)})
      
      output$boxplot <- renderPlot({ 
        boxplot(expr.matrix(), col= color(), ylab = "expression level", las = 2, 
                cex.axis = input$cex.axis_boxplot, 
                names = colnames(expr.matrix()),
                ylim = input$ylim,cex=0.1)
        mtext("Labels", side=1, line=7)
        title("Boxplot")
      })
      
      output$heatmap1 <- renderPlot({
        heatmap.2(expr.matrix.signif(), col=col_clu, scale="row",
                  main = "Heatmap result", key=TRUE, symkey=FALSE, 
                  density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=1, 
                  dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                  lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2),
                  RowSideColors=RowSideColors())
      },height = 1000, width = 600)
      
      output$heatmap2 <- renderPlot({
        heatmap.2(expr.matrix.signif.means(), col=col_clu, scale="row",
                  main = "Heatmap result", key=TRUE, symkey=FALSE, 
                  density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8, 
                  dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                  lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2),
                  RowSideColors=RowSideColors())
      },height = 1000, width = 600)
      
      
      output$MDS <- renderPlot({
        plotMDS(expr.matrix(),cex=input$cex.MDS,col=color())
      })
      
      output$PCA <- renderPlot({
        PCA(expr.matrix(),names2(),color(),namesColor(),input$cex.PCA)
      })
        
      output$volcanoplot <- renderPlot({
        volcanoPlot(topTable3()[[selectedComparison()]],input$logFCVolcano,input$pvalueVolcano)
      })
      
      output$MA <- renderPlot({
        MAPlot(topTable3()[[selectedComparison()]],input$logFCVolcano,input$pvalueVolcano)
      })
      
      output$heatDiagram <- renderPlot({
        heatmap.2(resultsContrast(), col=col_clu, scale="row",
                  main = "Heatmap result", key=TRUE, symkey=FALSE, 
                  density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8, 
                  dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                  lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2),
                  RowSideColors=RowSideColors())
      },height = 1000, width = 600)
      
      ################################################################################
      ###                              TABLES                                      ###
      ################################################################################
      output$resultsSummary <- renderDataTable({
        print(resultsSummary())
      })
      output$summaryClustering <- renderDataTable({
        print(summaryClust()[[(input$selectedClustn)]][[2]])
      })
      output$plotCorrectAttribution <- renderPlot({
        plotCorrectAttribution(summaryClust())
      })
      output$plotCorrectSubsample <- renderPlot({
        plotCorrectSubsample(summaryClust())
      })
      output$topTable <- renderDataTable({
        print(results()[[selectedComparison()]])
      })
      
      output$pathwayTable <- renderDataTable({
        print(toptableNodes1())
      })
      
      output$topGO <- renderDataTable({
        print(topGO())
      })
      output$topModuleGO <- renderDataTable({
        print(topModuleGO())
      })
      
      output$selectedGO <- renderDataTable({
        print(resultsOntology())
      })
      
      output$resultsModuleOntology <- renderDataTable({
        print(resultsModuleOntology())
      })
      output$modulesTable <- renderDataTable({
        print(modulesTable())
      })
      output$selectedModule <- renderDataTable({
        print(resultsModule()[[1]])
      })
      output$selectedModuleGO <- renderDataTable({
        print(ontologyModule())
      })
      output$selectedComparison <- renderPrint({
      })
      output$topReactions <- renderDataTable({
        print(selectedTopReactions())
      })
      output$toptableNodes <- renderDataTable({
        #print(toptableNodes1())
      })
      output$rawTopTable <- renderDataTable({
        print(topTable3()[[selectedComparison()]][match(input$selectedGenes,as.character(topTable3()[[selectedComparison()]][,typeID()])),])
      })
      
      ################################################################################
      ###                              DOWNLOAD HANDLER                            ###
      ################################################################################
      
      
      output$downloadExprMatrix <- downloadHandler(
        filename = function() { paste0("expr_matrix", '.rds') },
        content = function(file) {
          saveRDS(expr.matrix(), file)
        }
      )
      output$downloadSelectedGO <- downloadHandler(
        filename = function() { paste0("selectedGO", '.rds') },
        content = function(file) {
          saveRDS(resultsOntology(), file)
        }
      )
      
      output$downloadSummary <- downloadHandler(
        filename = function() { paste0("resultsSummary_",possibleComp(), '.rds') },
        content = function(file) {
          saveRDS(resultsSummary(), file)
        }
      )
      
      output$downloadResultsSelectedComparison <- downloadHandler(
        filename = function() { paste0("results_",selectedComparison(), '.rds') },
        content = function(file) {
          saveRDS(results()[[selectedComparison()]], file)
        }
      )
      output$downloadAllResults <- downloadHandler(
        filename = function() { paste0("results_", '.rds') },
        content = function(file) {
          saveRDS(results(), file)
        }
      )
      output$downloadNodeDataComparison <- downloadHandler(
        filename = function() { paste("node_rnaseq_",selectedComparison(), '.csv', sep='') },
        content = function(file) {
          write.csv(networkSelectedComparison()$nodeData, file)}
      )
      
      output$downloadEdgeDataComparison <- downloadHandler(
        filename = function() { paste("edge_rnaseq_",selectedComparison(), '.csv', sep='') },
        content = function(file) {
          write.csv(networkSelectedComparison()$edgeData, file)}
      )
      output$downloadNodeSelectedGO <- downloadHandler(
        filename = function() { paste("node_rnaseq_",selectedComparison(),input$selectedOntology, '.csv', sep='') },
        content = function(file) {
          write.csv(networkSelectedOntology()$nodeData, file)}
      )
      
      output$downloadEdgeSelectedGO <- downloadHandler(
        filename = function() { paste("edge_rnaseq_",selectedComparison(),input$selectedOntology, '.csv', sep='') },
        content = function(file) {
          write.csv(networkSelectedOntology()$edgeData, file)}
      )
      
      output$downloadSelectedGO <- downloadHandler(
        filename = function() { paste0("selectedGO", '.rds') },
        content = function(file) {
          saveRDS(resultsOntology(), file)}
      )
      
      #WILL NEED TO BE MODIFIED, or MAKE ANOTHER FUNCTION WITH RESULTS INSTEAD OF CONTRASTS; NOW IT TAKES EVERY GENES 
      output$downloadColorMatrix <- downloadHandler(
        filename = function() { paste0("colorMatrix", '.txt') },
        content = function(file) {
        write.table(
          addColorsMatrix(resultsContrast(),selectedComparison()),
          file,sep="\t\t",quote=FALSE,col.names=FALSE)}
      )
      
      
      output$downloadColorMatrixResults <- downloadHandler(
        filename = function() { paste0("colorMatrixComparison", '.txt') },
        content = function(file) {
          write.table(
            addColorsMatrixResults(results(),selectedComparison()),
            file,sep="\t\t",quote=FALSE,col.names=FALSE)}
      )
      
      output$downloadTopGO <- downloadHandler(
        filename = function() { paste0("topGO", '.rds') },
        content = function(file) {
          saveRDS(topGO(), file)}
      )
      output$downloadModulesTable <- downloadHandler(
        filename = function() { paste0("topGO", '.rds') },
        content = function(file) {
          saveRDS(modulesTable(), file)}
      )
      output$downloadtopModuleGO <- downloadHandler(
        filename = function() { paste0("topGO", '.rds') },
        content = function(file) {
          saveRDS(ontologyModule()[[1]], file)}
      )
      output$downloadResultsModuleOntology <- downloadHandler(
        filename = function() { paste0("topGO", '.rds') },
        content = function(file) {
          saveRDS(resultsModuleOntology(), file)}
      )
      observe({        
        updateSelectInput(session, "selectedComparison", choices = possibleComp())
      })
      observe({        
        updateSelectInput(session, "selectedModule", choices = as.numeric(head(modulesTable()[,"Module"])))
      })
      observe({ 
        if(input$submit1 < 2){
          updateSelectInput(session, "selectedVariables", choices = comparisons(), selected = comparisons())
        } else {
          x = input$selectedVariables
          updateSelectInput(session, "selectedVariables", choices = comparisons(), selected = x)
        }
      })
    })
    
    ################################################################################
    ###                              METHYLATION                                 ###
    ################################################################################
    
    #NOTES: 
    #1- meth.diff is the percentage of methylation difference. I don't know it counts it when one is 0...
    observeEvent(input$submit21, {
      
      #== INPUTS
      
      meth = reactive({readRDS(input$File5[,"datapath"])})
     
      if (!is.null(input$File6)){
        print("Results loading...")
        myDiff = reactive({readRDS(input$File6[,"datapath"])})
        print(names(myDiff()))
      } 
      
      topReactions = reactive({NULL})
      modules_table = reactive({NULL})
      expr.toBind = reactive({NULL})
      
      qvalue = eventReactive(input$submit22,{input$qvalue})
      methDiffNeg = eventReactive(input$submit22,{input$methDiffNeg})
      methDiffPos = eventReactive(input$submit22,{input$methDiffPos})
      #== Calculation of the beta-values matrix
      
      
      methy = reactive({getData(meth())})
      bvalues = reactive({calculateBvalues(methy(),attributes(meth())$sample.ids)})
      namesBvalues = reactive({get_names(colnames(bvalues()))[[1]]})
      names.uniqueBvalues = reactive({get_names(colnames(bvalues()))[[2]]})
      bvalues2 = reactive({conditionsChoice(input$selectedVariablesMethyl,bvalues(),namesBvalues())})
      design = reactive({design_contrasts(bvalues2())})
      color = reactive({color_groups(rownames(design()),colnames(design()))})
      
      if(input$submit21 < 2 ){
        print("names")
        possibleCompMethyl = reactive({possibleComparisons(namesBvalues())})
      } else {
        print("myDiff")
        possibleCompMethyl = eventReactive(input$submit21,{names(myDiff())})
      }
      
      #== Annotations
      
      #The files are on ratking
      transcripts_coord = matrix(unlist(strsplit(
        system('wget http://ratking.duckdns.org/~simon/rn5/rn5_annotation.txt -q -O -',intern=TRUE),
        "\t")),ncol=12,byrow=TRUE)
      transcripts_coord = transcripts_coord[,1:4]
      colnames(transcripts_coord) = c("chr","start","end","ensembl_transcript_id")
      
      transcripts_symbols = read.csv("annotations/databases/rat_symbols.csv")
      annotation1 = annotate_ensembl(as.character(transcripts_symbols[,"ensembl_transcript_id"]))
      cpg_diff_annotation = annotation1[[1]]
      go=annotation1[[2]]
      genes_annotation_unique = annotation1[[3]]
      typeID=annotation1[[4]]
      comparisonsMethyl = reactive({get_comparisons(namesBvalues())})
      
      observeEvent(input$submit23,{
        print("Calculating differential analysis...")
        myDiff = methylDifferentialAnalysisAll(meth(),possibleCompMethyl(),namesBvalues())
      })
      

      #== WGCNA modules input
      
      if (!is.null(input$File7)){
        bnet = reactive({readRDS(input$File7[,"datapath"])})
        expr.toBind = reactive({exprToBind(bnet(),bvalues2(),typeID())})
      } 
      
      #== Results
      allCpgs = reactive({bind_cpg_pos(myDiff(),0,0,1)})
      
      
      cpg = reactive({bind_cpg_pos(myDiff(),methDiffNeg(),methDiffPos(),qvalue())})
      
      
      #FAIRE FONCTION POUR OBTENIR list_cpg_cleared
      selectedCpg=reactive({queryCpg(allCpgs(),input$selectedCpgs,input$selectedComparisonMethyl,transcripts_coord,transcripts_symbols,typeID())})
      
      selectedGeneCpg = reactive({queryGene(allCpgs(),input$selectedCpgGenes,input$selectedComparisonMethyl,transcripts_coord,"ensembl_transcript_id")})
      selectedGeneCpgDiff = reactive({merge(selectedGeneCpg(),allCpgs()[[input$selectedComparisonMethyl]],by="cpg_pos")})
      cpg_annotation_list = reactive({cpg_annotation(cpg(),transcripts_coord,transcripts_symbols,"ensembl_transcript_id")})
      
      #cpg_annotation_list = reactive({cpg_annotation(cpg(),transcripts_symbols,"ensembl_transcript_id")})
      cpg_annotated = reactive({cpg_annotation_list()[[1]]})
      list_cpg_outside_gene = reactive({cpg_annotation_list()[[2]]})
      cpgByGene = reactive({findCpgByGene(cpg_annotated()[[input$selectedComparisonMethyl]],"symbol")})
      
      resultsOntologyMethyl = reactive({geneOntology(cpg_annotated(),go,"ensembl_transcript_id")})
      
      #== MethylDiff annotations
      
      #== Significative genes for heatmap
      
      bvalues.heatmap = reactive({bvalues.signif(bvalues2(),cpg_annotated(),namesBvalues(),names.uniqueBvalues(),"ensembl_transcript_id",input$symbol)})
      bvalues.heatmap.all = reactive({bvalues.heatmap()[[1]]})
      bvalues.heatmap.means = reactive({bvalues.heatmap()[[2]]})
      
      ################################################################################
      ###                              OUTPUTS                                     ###
      ################################################################################
      
      ################################################################################
      ###                              PLOTS                                       ###
      ################################################################################
      output$test2 <- renderText({
        print(names(summaryClust()))
      })
      output$test3 <- renderText({
        print(input$selectedComparisonMethyl)
      })
      #output$venn <- renderPlot({vennDiagram(results_contrast(),cex=input$venn.cex)})
      
      output$boxplot2 <- renderPlot({ 
        boxplot(bvalues2(), col= color(), ylab = "expression level", las = 2, 
                cex.axis = input$cex.axis_boxplot2, 
                names = colnames(bvalues2()),
                ylim = input$ylim2,cex=0.1)
        mtext("Labels", side=1, line=7)
        title("Boxplot")
      })
      
      output$heatmap21 <- renderPlot({
        heatmap.2(bvalues.heatmap.all(), col=col_clu, scale="row",
                  main = "Heatmap result", key=TRUE, symkey=FALSE, 
                  density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=1, 
                  dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                  lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2)
        )
      },height = 1000, width = 600)
      
      output$heatmap22 <- renderPlot({
        heatmap.2(bvalues.heatmap.means(), col=col_clu, scale="row",
                  main = "Heatmap result", key=TRUE, symkey=FALSE, 
                  density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8, 
                  dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                  lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2))
      })
      
      
      output$MDS2 <- renderPlot({
        plotMDS(bvalues2(),cex=input$cex.MDS,col=color())
      })
      output$PCA2 <- renderPlot({
        PCASamples(meth())
      })
      output$volcanoplot2 <- renderPlot({
        difference = myDiff[[1]]$meth.diff
        qvalue = myDiff[[1]]$qvalue
        plot(difference,-log(qvalue))
        abline(h= -log(0.35))
      })
      output$tree <- renderPlot({
        clusterSamples(meth(), dist = "correlation", method = "ward", plot = TRUE)
      })
      output$correlation <- renderPlot({
        getCorrelation(meth(), plot = T)
      })
        
      
      #output$MA <- renderPlot({
      #  difference = myDiff[[1]]$meth.diff
      
      #})
      
      output$heatDiagram <- renderPlot({
        hd=as.matrix(cpg_annotated[[1]]$meth.diff)
        colnames(hd) = names(cpg_annotated)
        rownames(hd) = cpg_annotated[[1]][,input$symbol]
        if (ncol(hd) < 2){
          heatmap.2(cbind(hd,hd), col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE, 
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8, 
                    dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2))
        }else{
          heatmap.2(hd, col=col_clu, scale="row",
                    main = "Heatmap result", key=TRUE, symkey=FALSE, 
                    density.info="none", trace="none", cexRow=input$cexRow_heatmap,cex=0.8, 
                    dendrogram=input$dendrogram, Colv = TRUE, Rowv = TRUE,
                    lhei=c(0.5,input$row_heatmap_height),lwid=c(0.5,2))
          
          
        }
        
      },height = 1000, width = 600)
      ################################################################################
      ###                              TABLES                                      ###
      ################################################################################
      
      output$topTableMethyl <- renderDataTable({
        print(cpg_annotated()[[input$selectedComparisonMethyl]])
      })
      output$cpgByGene <- renderDataTable({
        print(cpgByGene())
      })
      
      output$topTableMethylOutGene <- renderDataTable({
        print(list_cpg_outside_gene()[[input$selectedComparisonMethyl]])
      })
      output$topGOMethyl <- renderDataTable({
        print(cbind(rownames(resultsOntologyMethyl()[[input$selectedComparisonMethyl]][[2]]),resultsOntologyMethyl()[[input$selectedComparisonMethyl]][[2]]))
      })
      
      output$selectedGOMethyl <- renderDataTable({
        print(resultsOntologyMethyl())
      })
      
      output$rawTopCpgGenes <- renderDataTable({
        print(selectedGeneCpgDiff())
      })
      output$rawTopCpG <- renderDataTable({
        print(selectedCpg()[[1]][[1]])
      })
      
     
      
      #output$moduleTableMethyl <- renderDataTable({
      #  print(modules_tableMethyl())
      #})
      #output$selectedModuleMethyl <- renderDataTable({
      #  print(resultsModuleMethyl())
      #})
      
      #output$topReactionsMethyl <- renderDataTable({
      #  print(topReactionsMethyl()[[selectedComparisonMethyl()]])
      #})
      #output$toptableNodesMethyl <- renderDataTable({
      #  print(toptableNodes1Methyl())
      #})
      
      ################################################################################
      ###                              DOWNLOAD HANDLER                            ###
      ################################################################################
      
      
      #output$downloadSelectedGOMethyl <- downloadHandler(
      #  filename = function() { paste0("selectedGOMethyl", '.rds') },
      # content = function(file) {
      #saveRDS(resultsOntologyMethyl(), file)
      #}
      #)
      
      #output$downloadResultsMethyl <- downloadHandler(
      #  filename = function() { paste0("resultsMethyl", '.rds') },
      #  content = function(file) {
      #saveRDS(results()[[selectedComparisonMethyl()]], file)
      #  }
      #)
      
      #output$downloadNodeDataComparisonMethyl <- downloadHandler(
      #  filename = function() { paste("node_methyl_",selectedComparisonMethyl(), '.csv', sep='') },
      #  content = function(file) {
      #    #write.csv(networkSelectedComparisonMethyl()$nodeData, file)}
      #)
      
      #output$downloadEdgeDataComparisonMethyl <- downloadHandler(
      #  filename = function() { paste("edge_methyl_",selectedComparison(), '.csv', sep='') },
      #  content = function(file) {
      #    #write.csv(networkSelectedComparisonMethyl()$edgeData, file)}
      #)
      #output$downloadNodeSelectedGOMethyl <- downloadHandler(
      #  filename = function() { paste("node_methyl_",selectedComparison(),input$selectedOntology, '.csv', sep='') },
      #  content = function(file) {
      #    write.csv(networkSelectedOntologyMethyl()$nodeData, file)}
      #)
      
      #output$downloadEdgeSelectedGOMethyl <- downloadHandler(
      #  filename = function() { paste("edge_methyl_",selectedComparison(),input$selectedOntology, '.csv', sep='') },
      #  content = function(file) {
      #    #write.csv(networkSelectedOntologyMethyl()$edgeData, file)}
      #)
      
      #output$downloadSelectedGOMethyl <- downloadHandler(
      #  filename = function() { paste0("selectedGOMethyl", '.rds') },
      #  content = function(file) {
      #    #saveRDS(resultsOntologyMethyl(), file)}
      #)
      
      modules_tableMethyl = reactive({NULL})
      observe({        
        updateSelectInput(session, "selectedModuleMethyl", choices = as.numeric(head(modules_tableMethyl()[,"Module"],n=input$headModuleMethyl)))
      })
      observe({        
        updateSelectInput(session, "selectedComparisonMethyl", choices = possibleCompMethyl())
      })
      observe({ 
        if(input$submit21 < 2){
          updateSelectInput(session, "selectedVariablesMethyl", choices = comparisonsMethyl(), selected = comparisonsMethyl())
        } else {
          x = input$selectedVariablesMethyl
          updateSelectInput(session, "selectedVariablesMethyl", choices = comparisonsMethyl(), selected = x)
        }
      })
      
    })
    ################################################################################
    ###                    RNASEQ & METHYLATION INTEGRATION                      ###
    ################################################################################
    
    cpg_results_merged = reactive({merge(cpg_annotated()[[input$selectedComparison]],results()[[input$selectedComparisonMethyl]],by="ensembl_transcript_id")})
    output$mergeTable <- renderDataTable({
      print(cpg_results_merged())
    })
    
  })
