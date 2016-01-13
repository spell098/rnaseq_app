
library(shiny)
shinyUI(fluidPage(
  titlePanel("Linear model analysis"),
  
  tabsetPanel(
    tabPanel("RNA-seq",
             sidebarLayout(
               sidebarPanel(
                 helpText("Options"),    
                 fileInput("File1", "Updload RNAseq data file"),
                 selectInput("fileContent","File Content",choices=c("Expression Matrix","Expression Set"),selected="Expression Matrix"),
                 tags$textarea(id = "geo_id", placeholder = 'GEO ID', rows = 1, ""),
                 fileInput("File2", "Optional: Upload RNAseq module file"),
                 fileInput("File3", "Optional: Upload RNAseq reactions list"),
                 selectInput("specie","Specie",choices=c("hsapiens","rnorvegicus")),
                 selectInput("specieEnsembl","specieEnsembl",choices=c("rnorvegicus_gene_ensembl","hsapiens_gene_ensembl")),
                 selectInput("groupBy",label="Make group by",choices=c("groups","fuzzy"),selected="groups"),
                 numericInput("noCluster",label= "Number of clusters for fuzzy c-means clustering",value=3),
                 
                 selectInput("selectedComparison", label = "Select the comparison of interest", choices = NULL),
                 selectInput("selectedVariables",label= "Select the variables to compare",choices=NULL,multiple=TRUE),
                 
                 actionButton("submit1", "Run!"),
                 
                 numericInput("pvalue", label = "Select a p-value",value = 0.01,max=1,min=0),
                 numericInput("logFCneg", label = "Select the negative logFC limit", value = -1.3),
                 numericInput("logFCpos", label = "Select the positive logFC limit", value = 1.3),
                 selectInput("adjust",label="p-value correction method",selected="BH",choices=c("no","none","BH","BY","holm")),
                 actionButton("submit2", "Update restriction parameters")
                 
                 ,witdth=2),
               
               mainPanel(
                 helpText("test"),
                 textOutput("text1"),
                 tabsetPanel(
                   tabPanel("Quality Control",
                            tabsetPanel(
                              tabPanel("Boxplot",
                                       sliderInput("cex.axis_boxplot",label = "Caracters size on the axis",min=0,max=1,value=0.6),
                                       sliderInput("ylim",label = "Limits of the y-axis", min = -30, max = 30, value = c(-15, 15),step=0.1),
                                       plotOutput("boxplot"),
                                       help('Save data with new groups'),
                                       downloadButton("downloadExprMatrix","Download data")
                              ),
                              tabPanel("Clustering",
                                       numericInput("selectedClustn",label = "select number of clusters for summary",value=3),
                                       dataTableOutput("summaryClustering"),
                                       numericInput("pval_clust",label="select a p-value for significant clustering group",value=0.2),
                                       plotOutput("plotCorrectAttribution"),
                                       plotOutput("plotCorrectSubsample")
                              ),
                              tabPanel("PCA",
                                       sliderInput("cex.PCA", label = "Labels size of the PCA plot", min=0.5,max = 2, value = 1),
                                       plotOutput("PCA",height="500px"),
                                       sliderInput("cex.MDS", label = "Labels size of the MDS plot", min=0.5,max = 2, value = 1),
                                       plotOutput("MDS")
                                      
                              ),
                              tabPanel("Volcano & MA plots",
                                       numericInput("pvalueVolcano",label = "p-val for volcano plot",value = 0.05),
                                       numericInput("logFCVolcano",label = "Absolute logFC for volcano plot", value = 1),
                                       plotOutput("volcanoplot"),
                                       plotOutput("MA")                              
                                       
                              )
                            )
                            #N'arrive pas a afficher tous les graphique en loop
                            
                   ),
                   
                   tabPanel("Heatmaps",  
                            selectInput("dendrogram", label = "Dendrograms",choices = c("column","row","both","none"),selected = "none"),
                            selectInput("scale", label = "Scale", choices = c("column","row"), selected="row"),
                            sliderInput("cexRow_heatmap",label = "Caracters size on the x-axis",min=0,max=1.5,value=0.5),
                            sliderInput("row_heatmap_height",label = "Row height",min=0,max=6,value=2,step=0.05),
                            tabsetPanel(
                              tabPanel("Individual samples", 
                                       plotOutput("heatmap1",width="100%",height="1000px")                              
                              ),
                              tabPanel("Grouped samples",
                                       plotOutput("heatmap2",width="100%",height="1000px")                              
                              ),
                              tabPanel("Comparison heatmap",
                                       sliderInput("heatdiagram_cex",label = "Cex heatdiagram",min=0,max=3,value=0.6,step=0.05),
                                       plotOutput("heatDiagram",width="100%",height="1000px")
                              )
                            )
    
                   ) ,
                   tabPanel("Results",
                            tabsetPanel(
                              tabPanel("Results summary",
                                      help('Download results all comparisons'),
                                      downloadButton('downloadAllResults', 'Download'),
                                      dataTableOutput("resultsSummary"),
                                      downloadButton("downloadSummary","Download"),
                                      tags$textarea(id = "selectedGene", placeholder = 'gene of interest', rows = 8, ""),
                                      plotOutput("boxplot_element")
                              ),
                              tabPanel("Results table",
                                       dataTableOutput("topTable"),
                                       downloadButton('downloadResultsSelectedComparison', 'Download table'),
                                       downloadButton('downloadNodeDataComparison','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparison','Download edges file for cytoscape'),                              
                                       downloadButton("downloadColorMatrixResults","Download")
                              ),
                              tabPanel("Venn",
                                       sliderInput("lfc", label = "Select log Fold Change (logFC) limit (absolute)", value = c(1) ,max=5,min=0,step=0.05),
                                       sliderInput("venn.cex", label = "Venn Diagram caracters size", max=2,min=0.5,value=1),
                                       plotOutput("venn"),
                                       downloadButton("downloadVenn","Download")
                              ),
                              tabPanel("Top gene Ontologies",
                                       dataTableOutput("topGO"),
                                       downloadButton("downloadTopGO","Download"),
                                       tags$textarea(id = "selectedOntology", placeholder = 'ontologies of interest', rows = 8, ""),
                                       
                                       dataTableOutput("selectedGO"),
                                       downloadButton('downloadSelectedGO', 'Download'),
                                       downloadButton('downloadNodeSelectedGO','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeSelectedGO','Download edges file for cytoscape')
                              ),
                              tabPanel("Top modules",
                                       actionButton("submitMethylModule","Calculate Modules"),
                                       helpText("This action may take a while..."),
                                       dataTableOutput("modulesTable"),
                                       downloadButton("downloadModulesTable","Download"),
                                       dataTableOutput("selectedModule"),
                                       downloadButton("downloadSelectedModule","Download"),
                                       dataTableOutput("topModuleGO"),
                                       downloadButton("downloadtopModuleGO","Download"),
                                       dataTableOutput("resultsModuleOntology"),
                                       downloadButton("downloadResultsModuleOntology","Download")
                              ),
                              tabPanel("Raw data",
                                       #don't display the whole thing, it's useless; make the user enter is genes of interest manually, and display only that
                                       tags$textarea(id = "selectedGenes", placeholder = 'genes of interest', rows = 8, ""),
                                       dataTableOutput("rawTopTable"),
                                       downloadButton("downloadRawTopTable","Download"),
                                       plotOutput("boxplot_element_raw")                                       
                              )
               
                            )       
                            #verbatimTextOutput("results"),
    
                   ),
                   tabPanel("Network Analysis",
                            textOutput("confirmFindReactions"),
                            numericInput("threshold","Adjacency threshold for network",value = 0.7),
                            actionButton("findTopReactions","Find top reactions"),
                            dataTableOutput("topReactions"),
                            tags$textarea(id = "selectedReactions", placeholder = 'Reactions of interest', rows = 3, ""),
                            
                            #selectInput("selectedReactions", label = "Select pathways of interest", choices = NULL),
                            plotOutput("pathway",width="100%",height="800px"),
                            dataTableOutput("pathwayTable")
                   )
                 ))
             ) ),   
    tabPanel("Methylation",
             sidebarLayout(
               sidebarPanel(
                 helpText("Option"),    
                 fileInput("File5", "Updload metylation table file"),
                 fileInput("File6", "Updload metylation results file"),
                 fileInput("File7", "Optional: Upload methylation module file"),
                 fileInput("File8", "Optional: Upload methylation reactions list"),
                 
                 selectInput("specieEnsembl","specieEnsembl",choices=c("rnorvegicus_gene_ensembl","hsapiens_gene_ensembl")),
                 selectInput("symbol","Type of gene symbols",choices=c("rgd_symbol","hgnc_symbol")),
                 
                 selectInput("selectedComparisonMethyl", label = "Select the comparison of interest", choices = NULL),
                 selectInput("selectedVariablesMethyl",label= "Select the variables to compare",choices=NULL,multiple=TRUE),
                 
                 actionButton("submit21", "Run!"),
                 numericInput("methDiffNeg", label = "Select negative limit",value = -25),
                 numericInput("methDiffPos", label = "Select positive limit",value = 25),
                 
                 numericInput("qvalue", label = "Select a q-value",value = 0.01,max=1,min=0),
                 
                 actionButton("submit22", "Update restriction parameters")
                 
                 ,witdth=2),
               
               mainPanel(
                 textOutput("test2"),   
                 textOutput("test3"),

                 actionButton("submit23","Calculate differential analysis!"),
                 tabsetPanel(
                   tabPanel("Quality Control",
                            tabsetPanel(
                              tabPanel("Boxplot",
                                       sliderInput("cex.axis_boxplot2",label = "Caracters size on the axis",min=0,max=1,value=0.6),
                                       sliderInput("ylim2",label = "Limits of the y-axis", min = 0, max = 1, value = c(0, 1),step=0.1),
                                       plotOutput("boxplot2")
                              ),
                              tabPanel("MDS",
                                       sliderInput("cex.MDS", label = "Labels size of the MDS plot", min=0.5,max = 2, value = 1),
                                       plotOutput("MDS2"),
                                       plotOutput("PCA2")
                              ),
                              tabPanel("Volcanoplots",
                                       plotOutput("volcanoplot2")
                              ),
                              tabPanel("Correlation",
                                       plotOutput("correlation")
                              ),
                              tabPanel("Tree",
                                       plotOutput("tree")
                              )
                            )
                            #N'arrive pas a afficher tous les graphique en loop
                            
                   ),
                   
                   tabPanel("Heatmaps",  
                            selectInput("dendrogram2", label = "Dendrograms",choices = c("column","row","both","none"),selected = "none"),
                            selectInput("scale2", label = "Scale", choices = c("column","row"), selected="row"),
                            sliderInput("cexRow_heatmap2",label = "Caracters size on the x-axis",min=0,max=1.5,value=0.5),
                            sliderInput("row_heatmap_height2",label = "Row height",min=0,max=6,value=2,step=0.05),
                            tabsetPanel(
                              tabPanel("Individual samples", 
                                       plotOutput("heatmap21",width="100%",height="1000px")
                              ),
                              tabPanel("Grouped samples",
                                       plotOutput("heatmap22",width="100%",height="1000px")
                              ),
                              tabPanel("Comparison heatmap",
                                       sliderInput("heatdiagram_cex",label = "Cex heatdiagram",min=0,max=3,value=0.6,step=0.05),
                                       plotOutput("heatDiagram2",width="100%",height="1000px")
                              )
                            )                            
                   ) ,
                   tabPanel("Results",
                            tabsetPanel(
                              tabPanel("Results table",
                                       dataTableOutput("topTableMethyl"),
                                       downloadButton('downloadResultsMethyl', 'Download'),
                                       downloadButton('downloadNodeDataComparisonMethyl','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparisonMethyl','Download edges file for cytoscape'),
                                       dataTableOutput("topTableMethylOutGene"),
                                       downloadButton('downloadResultsMethylOut', 'Download'),
                                       downloadButton('downloadNodeDataComparisonMethylOutGene','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparisonMethylOutGene','Download edges file for cytoscape'),
                                       
                                       dataTableOutput("cpgByGene"),
                                       downloadButton('downloadResultsByGene', 'Download'),
                                       downloadButton('downloadNodeDataComparisonMethylByGene','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeDataComparisonMethylByGene','Download edges file for cytoscape')
                                       
                                       
                              ),
                              tabPanel("Venn",
                                       sliderInput("lfcMethyl", label = "Select log Fold Change (logFC) limit (absolute)", value = c(1) ,max=5,min=0,step=0.05),
                                       sliderInput("venn.cexMethyl", label = "Venn Diagram caracters size", max=2,min=0.5,value=1),
                                       plotOutput("vennMethyl") 
                              ),
                              tabPanel("Top gene Ontologies",
                                       dataTableOutput("topGOMethyl")
                              ),
                              tabPanel("Selected ontologies results",
                                       numericInput("headGOMethyl", label = "Number of gene ontologies to display", value = 10),
                                       tags$textarea(id = "selectedOntology", placeholder = 'ontologies of interest', rows = 8, ""),
                                       #selectInput("selectedOntology",label = "Select ontologies of interest",choices = NULL),
                                       numericInput("thresholdMethyl","Adjacency threshold for network",value = 0.7),
                                       
                                       dataTableOutput("selectedGOMethyl"),
                                       downloadButton('downloadSelectedGOMethyl', 'Download'),
                                       downloadButton('downloadNodeSelectedGOMethyl','Download nodes file for cytoscape'),
                                       downloadButton('downloadEdgeSelectedGOMethyl','Download edges file for cytoscape')
                              ),
                              tabPanel("Top modules",
                                       actionButton("submitMethylModule","Calculate Modules"),
                                       helpText("This action may take a while..."),
                                       dataTableOutput("moduleTableMethyl")
                                       
                              ),
                              tabPanel("Selected modules results",
                                       numericInput("headModuleMethyl", label = "Number of Modules to display", value = 20),
                                       selectInput("selectedModuleMethyl", label = "Select the module of interest", choices = NULL, multiple=TRUE),
                                       dataTableOutput("selectedModuleMethyl")
                                       
                              ),
                              tabPanel("Raw CpGs",
                                       tags$textarea(id = "selectedCpgs", placeholder = 'CpGs of interest', rows = 8, ""),
                                       dataTableOutput("rawTopCpG")
                              ),
                              tabPanel("Raw methylated genes",
                                       tags$textarea(id = "selectedCpgGenes", placeholder = 'Genes of interest', rows = 8, ""),
                                       dataTableOutput("rawTopCpgGenes")
                              )
                              
                            )       
                            #verbatimTextOutput("results"),
                            
                   ),
                   tabPanel("Network Analysis",
                            textOutput("confirmFindReactions2"),
                            actionButton("findTopReactions2","Find top reactions"),
                            dataTableOutput("topReactions2"),
                            tags$textarea(id = "selectedReactions2", placeholder = 'Reactions of interest', rows = 3, ""),
                            
                            #selectInput("selectedReactions", label = "Select pathways of interest", choices = NULL),
                            plotOutput("pathway2",width="100%",height="800px"),
                            dataTableOutput("toptableNodes2")           
                   )
                 ))
             )
    ),
    tabPanel("RNAseq / Methylation integration",
      tabsetPanel(
        
        tabPanel("CpG/RNA",
                 dataTableOutput("mergeTable")
                 ),
        tabPanel("Any CpG/RNA"
                 ),
        tabPanel("Network"
        )
      )
    )
  )))