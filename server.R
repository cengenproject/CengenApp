#### CenGenAPP 2018-2020
#### Report bugs to: gabrielsantperebaro@gmail.com


library(BiocManager)
require(shiny)
library(shinyjs)
library(shinythemes)
library(ggridges)
library(DT)
library(cowplot)
require(dplyr)
require(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(expss)
library(plotly)
library(xlsx)
#library(cairo)
#library(MAST)
options(repos = BiocManager::repositories())
#setwd("/Users/gabrielsantpere/Dropbox/Scripts/RshinyDeveloping_April2020")
#source("Functions.R")

utr <-
  c(
    "WBGene00023498",
    "WBGene00023497",
    "WBGene00004397",
    "WBGene00006843",
    "WBGene00004010",
    "WBGene00006789",
    "WBGene00001135",
    "WBGene00001079",
    "WBGene00001135",
    "WBGene00006783",
    "WBGene00000501",
    "WBGene00006788",
    "WBGene00001555",
    "WBGene00206533",
    "WBGene00011964",
    "WBGene00018172",
    "WBGene00016259",
    "WBGene00023407"
  )

server <- function(input, output) {
  ### Plot features ----
  observeEvent(input$PlotButton, {
    withProgress(message = "Generating Plots...",
                 value = 0,
                 expr = {
                   if (exists("id1")) {
                     removeNotification(id1)
                   }
                   
                   if (input$dataset == "Neurons only") {
                     all <- allNeurons
                   }
                   if (input$dataset == "All cell types") {
                     all <- allCells
                   }
                   Idents(object = all) <- "Neuron.type"
                   ranges <- reactiveValues(x = NULL, y = NULL)
                   observeEvent(input$plotPermanent_dblclick, {
                     brush <- input$plotPermanent_brush
                     if (!is.null(brush)) {
                       ranges$x <- c(brush$xmin, brush$xmax)
                       ranges$y <- c(brush$ymin, brush$ymax)
                     } else {
                       ranges$x <- NULL
                       ranges$y <- NULL
                     }
                   })
                   
                   
                   plots <- input$Plots
                   if (plots == 'Metadata') {
                     varToPlot <- input$featuresMetadata
                     featuresInput <- function() {
                       if (is.null(all)) {
                         return(NULL)
                       } else{
                         withProgress(message = "Generating Features Plot...", value = 0, {
                           DimPlot2(
                             object = all,
                             reduction = "umap",
                             group.by = 'Neuron.type',
                             label = TRUE,
                             pt.size = 0.7,
                             repel = TRUE,
                             l1 = ranges$x,
                             l2 = ranges$y
                           ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                 ylim = ranges$y,
                                                                                 expand = FALSE)
                         })
                       }
                     }
                     
                     if (input$featuresMetadata == 'Neuron.type' &&
                         input$dataset == "Neurons only") {
                       Type <- input$Cluster
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn1 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Neuron.type",
                               cells.highlight = names(all$Neuron.type[all$Neuron.type == Type]),
                               cols.highlight = "red",
                               pt.size = 0.7,
                               l1 = ranges$x,
                               l2 = ranges$y
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn1()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("NeuronType-", Type, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn1(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     if (input$featuresMetadata == 'Detection' &&
                         input$dataset == "Neurons only") {
                       Det <- input$Cluster3
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn3 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Detection" ,
                               cells.highlight = names(all$Detection[all$Detection == Det]),
                               cols.highlight = "red",
                               pt.size = 0.7
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn3()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("Detection-", Det, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn3(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     if (input$featuresMetadata == 'Experiment' &&
                         input$dataset == "Neurons only") {
                       Exp <- input$Cluster2
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn2 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Experiment",
                               cells.highlight = names(all$Experiment[all$Experiment == Exp]),
                               cols.highlight = "red",
                               pt.size = 0.7
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn2()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("Experiment-", Exp, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn2(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     if (input$featuresMetadata2 == 'Neuron.type' &&
                         input$dataset == "All cell types") {
                       Type <- input$Cluster.2
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn1 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Neuron.type",
                               cells.highlight = names(all$Neuron.type[all$Neuron.type == Type]),
                               cols.highlight = "red",
                               pt.size = 0.7,
                               l1 = ranges$x,
                               l2 = ranges$y
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn1()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("NeuronType-", Type, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn1(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     if (input$featuresMetadata2 == 'Detection' &&
                         input$dataset == "All cell types") {
                       Det <- input$Cluster3.2
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn3 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Detection" ,
                               cells.highlight = names(all$Detection[all$Detection == Det]),
                               cols.highlight = "red",
                               pt.size = 0.7
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn3()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("Detection-", Det, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn3(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     if (input$featuresMetadata2 == 'Experiment' &&
                         input$dataset == "All cell types") {
                       Exp <- input$Cluster2.2
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn2 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Experiment",
                               cells.highlight = names(all$Experiment[all$Experiment == Exp]),
                               cols.highlight = "red",
                               pt.size = 0.7
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn2()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("Experiment-", Exp, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn2(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     if (input$featuresMetadata2 == 'Tissue.type' &&
                         input$dataset == "All cell types") {
                       output$SinglePlot <-
                         renderPlot({
                           featuresInput()
                         }, height = 400, width = 600)
                       fn4 <-
                         function() {
                           withProgress(message = "Generating Features Plot...", value = 0, {
                             DimPlot2(
                               object = all,
                               reduction = "umap",
                               group.by = "Tissue.type",
                               cells.highlight = names(all$Tissue.type[all$Tissue.type == input$Cluster4]),
                               cols.highlight = "red",
                               pt.size = 0.7
                             ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                                   ylim = ranges$y,
                                                                                   expand = FALSE)
                           })
                         }
                       output$High <-
                         renderPlot({
                           fn4()
                         }, height = 400, width = 600, execOnResize = FALSE)
                       output$downloadHigh <-
                         downloadHandler(
                           filename = function() {
                             filename = paste("Tissue-", input$Cluster4, ".png", sep = "")
                           },
                           content = function(file) {
                             ggsave(
                               fn4(),
                               file = file,
                               height = 100,
                               width = 150,
                               units = "mm" ,
                               limitsize = FALSE,
                               device = "png"
                             )
                           }
                         )
                     }
                     
                     output$downloadClust <-
                       downloadHandler(
                         filename = function() {
                           filename = paste("Clusters.png", sep = "")
                         },
                         content = function(file) {
                           ggsave(
                             featuresInput(),
                             file = file,
                             height = 200,
                             width = 300,
                             units = "mm" ,
                             limitsize = FALSE,
                             device = "png"
                           )
                         }
                       )
                   }
                   
                   ### Plot genes individual ----
                   if (plots == 'Individual genes') {
                     g <- input$GeneName
                     print(g)
                     
                     if (g %in% gene_list$gene_name) {
                       varToPlot <- filter(gene_list, gene_name == g)$gene_id
                     }
                     
                     if (g %in% gene_list$gene_id) {
                       varToPlot <- g
                     }
                     
                     if (g %in% gene_list$seqnames) {
                       varToPlot <- filter(gene_list, seqnames == g)$gene_id
                     }
                     
                     output$geneUTR <-
                       isolate(renderText({
                         shiny::validate(need(
                           !varToPlot %in% utr,
                           message = paste0(
                             "WARNING: ",
                             g,
                             " expression is unreliable as it has been overexpressed to generate transgenic strains."
                           )
                         ))
                       }))
                     
                     
                     f <- function() {
                       DimPlot2(
                         object = all,
                         reduction = "umap",
                         group.by = "Neuron.type",
                         label = TRUE,
                         repel = TRUE,
                         l1 = ranges$x,
                         l2 = ranges$y
                       ) + theme(legend.position = "none") + coord_cartesian(xlim = ranges$x,
                                                                             ylim = ranges$y,
                                                                             expand = TRUE)
                     }
                     f2 <- function() {
                       RidgePlot(
                         object = all,
                         features = as.character(varToPlot),
                         ncol = 1,
                         group.by = "Neuron.type",
                         sort = "decreasing"
                       ) + theme(
                         legend.position = "none",
                         panel.grid = element_blank(),
                         panel.background = element_blank()
                       )
                     }
                     f3 <-
                       function() {
                         FeaturePlot3(
                           object = all,
                           features = as.character(varToPlot),
                           reduction = "umap",
                           pt.size = 0.7,
                           combine = TRUE,
                           cols = c("beige", "darkred"),
                           coord.fixed = FALSE,
                           ranges = ranges
                         )
                       }
                     
                     output$SinglePlot2 <- renderPlot({
                       shiny::validate(
                         need(as.character(g) %in% gene_list$gene_name | as.character(g) %in% gene_list$gene_id | as.character(g) %in% gene_list$seqnames,  message = "Gene name is not present in this dataset or wrong name")
                       )
                       f()
                     }, height = 400, width = 600, execOnResize = FALSE)
                     output$FeaturePlot <-  renderPlot({
                       shiny::validate(
                         need(as.character(g) %in% gene_list$gene_name | as.character(g) %in% gene_list$gene_id | as.character(g) %in% gene_list$seqnames,  message = "Gene name is not present in this dataset or wrong name")
                       )
                       f3()
                     }, height = 400, width = 600, execOnResize = FALSE)
                     output$Violin <- renderPlot({
                       shiny::validate(
                         need(as.character(g) %in% gene_list$gene_name | as.character(g) %in% gene_list$gene_id | as.character(g) %in% gene_list$seqnames,  message = "Gene name is not present in this dataset or wrong name")
                       )
                       f2()
                     }, height = 1300, width = 600)
                     
                     
                     output$downloadExp <-
                       downloadHandler(
                         filename = function() {
                           filename = paste("Dim-", varToPlot, ".png", sep = "")
                         },
                         content = function(file) {
                           ggsave(
                             f3(),
                             file = file,
                             height = 100,
                             width = 150,
                             units = "mm" ,
                             limitsize = FALSE,
                             device = "png"
                           )
                         }
                       )
                     output$downloadViolin <-
                       downloadHandler(
                         filename = function() {
                           filename = paste("Violin-", varToPlot, ".png", sep = "")
                         },
                         content = function(file) {
                           ggsave(
                             f2(),
                             file = file,
                             height = 1300,
                             width = 600,
                             units = "mm" ,
                             limitsize = FALSE,
                             device = "png"
                           )
                         }
                       )
                   }
                   
                   ### Plot genes colocoalization ----
                   
                   if (plots == 'Colocalization') {
                     g1 <- input$ColocalizationGenes1
                     g2 <- input$ColocalizationGenes2
                     #gene1 <- filter(gene_list, gene_name == g1)$gene_id
                     #gene2 <- filter(gene_list, gene_name == g2)$gene_id
                     
                     if (g1 %in% gene_list$gene_name) {
                       varToPlot1 <- filter(gene_list, gene_name == g1)$gene_id
                     }
                     
                     if (g1 %in% gene_list$gene_id) {
                       varToPlot1 <- g1
                     }
                     
                     if (g1 %in% gene_list$seqnames) {
                       varToPlot1 <- filter(gene_list, seqnames == g1)$gene_id
                     }
                     
                     if (g2 %in% gene_list$gene_name) {
                       varToPlot2 <- filter(gene_list, gene_name == g2)$gene_id
                     }
                     
                     if (g2 %in% gene_list$gene_id) {
                       varToPlot2 <- g2
                     }
                     
                     if (g2 %in% gene_list$seqnames) {
                       varToPlot2 <- filter(gene_list, seqnames == g2)$gene_id
                     }
                     
                     
                     featuresInput <- function() {
                       if (is.null(all)) {
                         return(NULL)
                       } else{
                         varBlend <- input$blend
                         withProgress(message = "Generating Features Plot...", value =
                                        0, {
                                          FeaturePlot2(
                                            object = all,
                                            features = c(
                                              #as.character(filter(gene_list, gene_name == g1)$gene_id),
                                              #as.character(filter(gene_list, gene_name == g2)$gene_id)
                                              varToPlot1, varToPlot2
                                            ),
                                            reduction = "umap",
                                            pt.size = 0.7,
                                            combine = TRUE,
                                            cols = c("beige", "darkred"),
                                            blend = TRUE,
                                            blend.threshold = varBlend,
                                            coord.fixed = TRUE,
                                            ranges = ranges
                                          )
                                        })
                       }
                     }
                     
                     output$SinglePlotDouble <- isolate(renderPlot({
                       shiny::validate(
                         need(as.character(g1) %in% gene_list$gene_name | as.character(g1) %in% gene_list$gene_id | as.character(g1) %in% gene_list$seqnames,  message = "Gene name 1 is not present in this dataset or wrong name")
                       )
                       shiny::validate(
                         need(as.character(g2) %in% gene_list$gene_name | as.character(g2) %in% gene_list$gene_id | as.character(g2) %in% gene_list$seqnames,  message = "Gene name 2 is not present in this dataset or wrong name")
                       )
                       featuresInput()
                     }, height = 400, width = 1000))
                     output$gene1utr <-
                       renderText({
                         shiny::validate(need(
                           !varToPlot1 %in% utr ,
                           message = paste0(
                             "WARNING: ",
                             g1,
                             " expression is unreliable as it has been overexpressed to generate transgenic strains."
                           )
                         ))
                       })
                     output$gene2utr <-
                       renderText({
                         shiny::validate(need(
                           !varToPlot2 %in% utr ,
                           message = paste0(
                             "WARNING: ",
                             g2,
                             " expression is unreliable as it has been overexpressed to generate transgenic strains."
                           )
                         ))
                       })
                     output$downloadCol <-
                       downloadHandler(
                         filename = function() {
                           filename = paste("Colocalization-", g1, "_", g2, ".png", sep = "")
                         },
                         content = function(file) {
                           ggsave(
                             featuresInput(),
                             file = file,
                             height = 200,
                             width = 500,
                             units = "mm" ,
                             limitsize = FALSE,
                             device = "png"
                           )
                         }
                       )
                     output$SinglePlotPermanent <-
                       renderPlot(
                         DimPlot2(
                           object = all,
                           reduction = "umap",
                           group.by = "Neuron.type",
                           label = TRUE,
                           repel = TRUE,
                           l1 = ranges$x,
                           l2 = ranges$y
                         ) + theme(legend.position = "none") + coord_cartesian(
                           xlim = ranges$x,
                           ylim = ranges$y,
                           expand = FALSE
                         ),
                         height = 400,
                         width = 600,
                         execOnResize = FALSE
                       )
                     
                   }
                 })
    
  })
  
  ### Tables of markers ----
  markers$p_val <-
    formatC(markers$p_val, format = "e", digits = 3) %>% gsub(" ", "", .)
  markers$p_val_adj <-
    formatC(markers$p_val_adj, format = "e", digits = 3) %>% gsub(" ", "", .)
  markers$avg_log2FC <-
    formatC(markers$avg_log2FC, digits = 3) %>% gsub(" ", "", .)
  markersAllcells$p_val <-
    formatC(markersAllcells$p_val,
            format = "e",
            digits = 3) %>% gsub(" ", "", .)
  markersAllcells$p_val_adj <-
    formatC(markersAllcells$p_val_adj,
            format = "e",
            digits = 3) %>% gsub(" ", "", .)
  markersAllcells$avg_log2FC <-
    formatC(markersAllcells$avg_log2FC, digits = 3) %>% gsub(" ", "", .)
  
  observeEvent(input$Markers, {
    print(input$top2)
    output$MarkTable <- DT::renderDataTable({
      DT::datatable(
        filter(markers, cluster == input$Markers, avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% head(as.numeric(input$top)) ,
        style = 'jQueryUI',
        class = 'cell-border stripe',
        rownames = FALSE
      ) %>% formatStyle(c(1:8), color = "black")
    })
    t1 <-
      filter(markers, cluster == input$Markers, avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC))
    output$downloadMarkers <-
      downloadHandler(
        filename = function() {
          paste("CellMarkers_NeuronsOnly-", input$Markers, ".csv", sep = "")
        },
        content = function(file) {
          write.csv(t1, file, sep = "\t")
        }
      )
  })
  observeEvent(input$Markers2, {
    output$MarkTable2 <- DT::renderDataTable({
      DT::datatable(
        filter(markersAllcells, cluster == input$Markers2, avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% head(as.numeric(input$top2)),
        style = 'jQueryUI',
        class = 'cell-border stripe',
        rownames = FALSE
      ) %>% formatStyle(c(1:8), color = "black")
    })
    t2 <-
      filter(markersAllcells, cluster == input$Markers2, avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC))
    output$downloadMarkers2 <-
      downloadHandler(
        filename = function() {
          paste("CellMarkers_AllCells-", input$Markers2, ".csv", sep = "")
        },
        content = function(file) {
          write.csv(t2, file, sep = "\t")
        }
      )
  })
  
 
  observeEvent(input$DEXButton, {
    print(input$batch1)
    print(input$batch2)
    Idents(object = allCells) <- "Neuron.type"
    b1 <- unlist(strsplit(input$batch1, split = ","))
    b1 <- gsub(" ", "", as.character(b1))
    b2 <- unlist(strsplit(input$batch2, split = ","))
    b2 <- gsub(" ", "", as.character(b2))
    
    if (mean(c(b1, b2) %in% allCells$Neuron.type) == 1 | b2 == "ALL" | b2 == "NEURONS") {
      output$text2 <- renderText({
        ""
      })
      
      if (!b2 %in% c("ALL", "NEURONS")){
      
        tableDEX <-
          FindMarkers(
            allCells,
            #ident.1 = unlist(strsplit(input$batch1, split = ",| ")),
            #ident.2 = unlist(strsplit(input$batch2, split = ",| ")),
            ident.1 = b1,
            ident.2 = b2,
            test.use = input$test
          )
        tableDEX$gene <- rownames(tableDEX)
      
      }
      
      if (b2 == "ALL"){
        print(b2)
         tableDEX <-
          FindMarkers(
            allCells,
            #ident.1 = unlist(strsplit(input$batch1, split = ",| ")),
            ident.1 = b1,
            test.use = input$test
          )
        
         tableDEX$gene <- rownames(tableDEX)
      }
      
      if (b2 == "NEURONS"){
        print(b2)
        tableDEX <-
          FindMarkers(
            allCells,
            #ident.1 = unlist(strsplit(input$batch1, split = ",| ")),
            ident.1 = b1,
            ident.2 = filter(allCells@meta.data, !Neuron.type %in% b1, Tissue.type == "Neuron")$Neuron.type %>% unique,
            test.use = input$test
          )
        
        tableDEX$gene <- rownames(tableDEX)
      }
      
      
        tableDEX <-
        merge(
          tableDEX,
          gene_list,
          by.x = "gene",
          by.y = "gene_id",
          all.x = TRUE
        )
      
      if (input$test != "roc") {
        tableDEX <- tableDEX %>% arrange(p_val_adj)
        tableDEX$p_val <-
          as.numeric(formatC(tableDEX$p_val, format = "e", digits = 3) %>% gsub(" ", "", .))
        tableDEX$p_val_adj <-
          as.numeric(formatC(
            tableDEX$p_val_adj,
            format = "e",
            digits = 3
          ) %>% gsub(" ", "", .))
      }
      if (input$test == "roc") {
        tableDEX <- tableDEX %>% arrange(desc(avg_log2FC))
      }    
        tableDEX$avg_log2FC <-
          as.numeric(formatC(tableDEX$avg_log2FC, digits = 3) %>% gsub(" ", "", .))
  
        
        
      if (nrow(tableDEX) > 0) {
        output$MarkTable_Batch <- DT::renderDataTable({
          DT::datatable(
            tableDEX %>% head(input$topM2),
            options = list(pageLength = input$topM2),
            style = 'jQueryUI',
            class = 'cell-border stripe',
            rownames = FALSE
          ) %>% formatStyle(c(1:9), color = "black")
        })
        output$downloadDEX <-
          downloadHandler(
            filename = function() {
              paste(
                "DEXGens-",
                paste(b1, collapse = ","),
                "-",
                paste(b2, collapse = ","),
                ".csv",
                sep = ""
              )
            },
            content = function(file) {
              write.csv(tableDEX, file, sep = "\t")
            }
          )
      } else {
        output$MarkTable_Batch <-
          "No features pass logfc.threshold threshold"
      }
    } else {
      output$MarkTable_Batch <- DT::renderDataTable({
        "NULL"
      })
      output$text2 <-
        renderText({
          "One or more cell types introduced are not correct"
        })
    }
  }, ignoreNULL = TRUE)
  
  ### React to thresholds ----
  
  observeEvent(input$TCell, {
    withProgress(message = "Obtaining information...", value = 0, {
      
      if( input$Tcell_cut == "All Cells Unfiltered" ) { th = L4.all.TPM.raw_th } else { th = ths }
      
      output$Error1 <-
        isolate(renderText({
          shiny::validate(need(
            input$Tcell_name %in% colnames(th),
            message = paste0(
              "WARNING: If you want to query non-neuronal cells select All Cells Unfiltered"
            )
          ))
        }))
      
      if (input$Tcell_name %in% colnames(th)) {
        t4 <- dplyr::filter(th, threshold == input$Tcell_cut) %>% dplyr::select(gene_name, input$Tcell_name)
        t4d <- dplyr::filter(th, threshold == input$Tcell_cut) %>% dplyr::select(gene_name, X, input$Tcell_name)
        t4 <- t4[rev(order(t4[, 2])), ]
        t4d <- t4d[rev(order(t4d[, 3])), ]
        t4[, 2] <- as.numeric(formatC(t4[, 2], digits = 3, format = "f") %>% gsub(" ", "", .))
        colnames(t4) <- c("Gene name", "Expression level")
        t4d[, 3] <- as.numeric(formatC(t4d[, 3], digits = 3, format = "f") %>% gsub(" ", "", .))
        colnames(t4d) <- c("Gene name", "Gene ID", "Expression level")
      
       output$Tcell_name_table <-
         DT::renderDataTable({
           DT::datatable(
            t4,
            options = list(pageLength = 10, autoWidth = TRUE),
            rownames = FALSE,
            style = 'jQueryUI',
            class = 'cell-border stripe'
          ) %>% formatStyle(c(1:2), color = "black")
        })
      
      output$get_download_gene <- renderUI({   
        req(input$TCell)
        downloadButton('downloadGene', "Download table")
        })
      
      output$downloadGene <-
        downloadHandler(
          filename = function() {
            paste(
              "GenesExpressed_in_",
              input$Tcell_name,
              "-thrs",
              input$Tcell_cut,
              ".csv",
              sep = ""
            )
          },
          content = function(file) {
            write.csv(t4d , file, dec = ".", sep = "\t")
          }
        )
    
      } else { output$Tcell_name_table  <- NULL }
      
  })
  })
  
  observeEvent(input$TGene, {
    
    g <- input$Tgene_name
    print(g)
    
    var <- ""
    
    if (g %in% gene_list$gene_name) {
      var <- g
    }
    
    if (g %in% gene_list$gene_id) {
      var <- dplyr::filter(gene_list, gene_id == g)$gene_name
    }
    
    if (g %in% gene_list$seqnames) {
      var <- dplyr::filter(gene_list, seqnames == g)$gene_name
    }
    
    if( input$Tgene_cut == "All Cells Unfiltered" ) { 
      th = L4.all.TPM.raw_th
      columns <- c(1:169)
    } else { 
      th = ths
      columns <- c(2:129)
    }
    withProgress(message = "Obtaining information...", value = 0, {
      if (var %in% dplyr::filter(th, threshold == input$Tgene_cut)$gene_name) {
        t3 <-
          sort(
            filter(
              th,
              gene_name == var,
              threshold == input$Tgene_cut
            )[, columns],
            decreasing = TRUE
          )
        t3 <- data.frame(CellType = names(t3), expression = as.numeric(t3))
        t3 <- dplyr::filter(t3, expression > 0)
        if (nrow(t3) > 0) {
          t3[, 2] <-
            as.numeric(formatC(t3[, 2], digits = 3, format = "f") %>% gsub(" ", "", .))
        }
        colnames(t3) <- c("Cell type", "Expression level")
        output$text1 <-
          isolate(renderText({
            shiny::validate(need(
              !dplyr::filter(gene_list, gene_name == var)$gene_id %in% utr,
              message = paste0(
                "WARNING: ",
                input$Tgene_name,
                " expression is unreliable as it has been overexpressed to generate transgenic strains."
              )
            ))
          }))
        #output$text1 <- renderText({""})
      }
      else {
        print("NO")
        t3 <- NULL
        output$text1 <-
          renderText({
            "Gene is not expressed at any threshold or does not exist"
          })
      }
      
      output$Tgene_name_table <-
        DT::renderDataTable({
          DT::datatable(
            t3,
            options = list(pageLength = 10, autoWidth = TRUE),
            rownames = FALSE,
            style = 'jQueryUI',
            class = 'cell-border stripe'
          ) %>% formatStyle(c(1:2), color = "black")
        })
      
      output$get_download_cell <- renderUI({   
        req(input$Tgene_name)
        downloadButton('downloadCell', "Download table")
      })
      
      output$downloadCell <-
        downloadHandler(
          filename = function() {
            paste(
              "GenesExpressing-",
              input$Tgene_name,
              "-thrs",
              input$Tgene_cut,
              ".csv",
              sep = ""
            )
          },
          content = function(file) {
            write.csv(t3, file, dec = ".", sep = "\t")
          }
        )
    })
  })
  
  observeEvent( list(input$Tgene_name_batch, input$Tgene_cut_batch) , {
    
    gns1 <- strsplit(as.character(input$Tgene_name_batch), "\n| |\\,|\t")
    gns1 <- as.data.frame(gns1)[,1]
    gns <- unique(c(gns1, filter(gene_list, gene_id %in% gns1 | seqnames %in% gns1)$gene_name))
    
    if( input$Tgene_cut_batch == "All Cells Unfiltered" ) { 
      th = L4.all.TPM.raw_th
      columns <- c(171, 170, 172, c(1:169) )
    } else { 
      th = ths
      columns <- c(1, 131, 130, c(2:129))    
    }
    
    if (length(which(gns %in% unique(th$gene_name))) > 0) {
      tb <-
        dplyr::filter(th,
               gene_name %in% gns,
               threshold == input$Tgene_cut_batch)[, columns]
      output$textb <- renderText({""})
      head(tb)
      req(input$Tgene_name_batch)
      
      output$TGeneBatch <-
        downloadHandler(
          filename = function() {
            paste(
              "GenesExpressing-BATCH",
              "-thrs",
              input$Tgene_cut_batch,
              ".xlsx",
              sep = ""
            )
          },
          content = function(file) {
            write.xlsx(tb, file, sheetName = "Sheet1", 
                       col.names = TRUE, row.names = TRUE, append = FALSE)
          }
        )
      
    } else {
      print("NO")
      tb <- NULL
      output$textb <-
        renderText({
          "No gene found with this name in the dataset"
        })
    }
  })
  
  ### Percentages ----
  
  observeEvent(input$Filter, {
    s1 <- unlist(strsplit(as.character(input$String1), split = ","))
    s1 <- gsub(" ", "", as.character(s1))
    s2 <- unlist(strsplit(as.character(input$String2), split = ","))
    s2 <- gsub(" ", "", as.character(s2))
    expressed <- input$Expressed
    not_expressed <- input$NonExpressed
    
    print(expressed)
    yes <-
      filter(pcttable, id %in% s1, pct.exp >= as.numeric(expressed))
    yes.genes <- names(which(table(yes$gene_name) == length(s1)))
    no <-
      filter(pcttable, id %in% s2, pct.exp < as.numeric(not_expressed))
    no.genes <- names(which(table(no$gene_name) == length(s2)))
    tt3 <- intersect(yes.genes, no.genes)
    tt3t <- unique(filter(pcttable, gene_name %in% tt3)[, c(1, 5)])
    tt1 <- filter(yes, gene_name %in% tt3)
    tt2 <- filter(no, gene_name %in% tt3)
    
    
    output$YesExpressed <-
      DT::renderDataTable({
        DT::datatable(
          yes,
          options = list(pageLength = 10, autoWidth = TRUE),
          rownames = FALSE,
          style = 'jQueryUI',
          class = 'cell-border stripe'
        ) %>% formatStyle(c(1:6), color = "black")
      })
    output$NoExpressed <-
      DT::renderDataTable({
        DT::datatable(
          no,
          options = list(pageLength = 10, autoWidth = TRUE),
          rownames = FALSE,
          style = 'jQueryUI',
          class = 'cell-border stripe'
        ) %>% formatStyle(c(1:6), color = "black")
      })
    output$Result <-
      DT::renderDataTable({
        DT::datatable(
          as.data.frame(tt3t),
          options = list(pageLength = 10, autoWidth = TRUE),
          rownames = FALSE,
          style = 'jQueryUI',
          class = 'cell-border stripe'
        ) %>% formatStyle(c(1:2), color = "black")
      })
    
    output$downloadQuery <-
      downloadHandler(
        filename = function() {
          paste(
            "GenesProportion-",
            paste(s1, collapse = ","),
            "-thrs",
            paste(s2, collapse = ","),
            ".csv",
            sep = ""
          )
        },
        content = function(file) {
          write.csv(rbind(tt1, tt2), file, sep = "\t")
        }
      )
    
  })


#############################################################################
########### Heatmaps ########################################################
#############################################################################


observeEvent(input$PlotHeatmap, {
  ds <- input$dataset_heatmap
  #ss <- unlist(strsplit(as.character(input$genelist), split = ","))
  #ss <- gsub(" ", "", as.character(ss))
  ss <- strsplit(as.character(input$genelist), "\n| |\\,|\t")
  ss <- as.data.frame(ss)[,1]
  ss <- unique(c(ss, filter(gene_list, gene_id %in% ss | seqnames %in% ss)$gene_name))
 
  mis <- ss[ss %in% c(gene_list$gene_id, gene_list$gene_name, gene_list$seqnames) & !ss %in% med.scaled.long$gene_name & ss %in% gene_list$gene_name]
  mis_all <- ss[ss %in% c(gene_list$gene_id, gene_list$gene_name, gene_list$seqnames) & !ss %in% L4.TPM.raw.scaled.long$gene_name & ss %in% gene_list$gene_name]
  
  if(ds=="Neurons only"){
    L4.TPM=L4.TPM.medium
    heatmapdata=med.scaled.long
    cc = colnames(ths)[-c(1,130,131)]
    missing = mis} else {
    L4.TPM=as(L4.all.TPM.raw,"dgCMatrix")
    heatmapdata=L4.TPM.raw.scaled.long
    cc=colnames(L4.all.TPM.raw)
    missing = mis_all
    }
  
  head(heatmapdata)
  print(ds)
  
  flp.neuron.scaled <- heatmapdata[which(heatmapdata$gene_name %in% ss),]
  flp.ids <- as.character(vlookup(unique(flp.neuron.scaled$gene_name), gene_list, result_column = 1, lookup_column = 2))
  flp.expr <- L4.TPM[flp.ids,, drop=FALSE]
  if ( nrow(flp.expr) >1 ) {
    flp.neuron.order <- pheatmap(flp.expr, scale = "row")
    flp.neuron.order <- flp.neuron.order[["tree_row"]]$order
    flp.neuron.order <- flp.ids[flp.neuron.order]
  } else {
    flp.neuron.order <- flp.ids[1]
  }
  flp.neuron.order <- as.character(vlookup(flp.neuron.order, gene_list))
  flp.neuron.scaled$gene_name <- factor(flp.neuron.scaled$gene_name, levels = c(rev(flp.neuron.order), missing))
  #flp.neuron.scaled$gene_name <- fct_rev(flp.neuron.scaled$gene_name)
  
  for( i in missing ){
    dff <- data.frame(gene_name=i, cell.type= cc, scaled.expr=0, prop=0, Modality="NA")
    if(ds!="Neurons only"){colnames(dff)[5]<-"tissue"}
    flp.neuron.scaled <- rbind(flp.neuron.scaled, dff)
  }
  
  if (ds!="Neurons only"){
    
   if(nrow(flp.expr) >1){
    g <- ggplot(flp.neuron.scaled, aes(y = gene_name, x = cell.type)) + 
      geom_point(aes(color = scaled.expr, size = prop)) +
      theme(panel.background = element_blank(), axis.text.y.left = element_text(size = 12, color = "black"),
            legend.key = element_blank(),
            legend.text = element_text(color = "black", size = 14),
            legend.title = element_text(color = "black", size = 16),
            axis.title = element_text(size = 16, color = "black")) + 
      scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) + 
      scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) + 
      labs(y = "Gene", x = "Tissue") + theme(panel.grid = element_line(size = 0.2, color = "grey85")) +
      facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") + 
      theme(strip.placement = "outside", 
            strip.background.x = element_blank(),
            axis.text.x.bottom = element_blank(),
            strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) 
   } else {
     
     g <- ggplot(flp.neuron.scaled, aes(y = gene_name, x = cell.type)) + 
       geom_point(aes(color = scaled.expr, size = prop)) +
       theme(panel.background = element_blank(), axis.text.y.left = element_text(size = 12, color = "black"),
             legend.key = element_blank(),
             legend.text = element_text(color = "black", size = 14),
             legend.title = element_text(color = "black", size = 16),
             axis.title = element_text(size = 16, color = "black")) + 
       scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) + 
       scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) + 
       labs(y = "Gene", x = "Tissue") + theme(panel.grid = element_line(size = 0.2, color = "grey85"))
   }
    
    pg <- ggplotGrob(g)
    
    for(i in which(grepl("strip-b", pg$layout$name))){
      pg$grobs[[i]]$layout$clip <- "off"
    }
    
  fnh <-
    function() {
      withProgress(message = "Generating heatmap Plot...", value = 0, {
        grid::grid.draw(pg)      })
    }
  } else {
  
    g<-ggplot(flp.neuron.scaled, aes(y =gene_name, x = cell.type)) +
      geom_point(aes(color = scaled.expr, size = prop)) +
      #geom_point(data = NULL, aes(y = vec), pch = NA) +
      theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 7, color = "black", hjust = 1),
            panel.background = element_blank(), axis.text.y.left = element_text(size = 7, color = "black")) +
      scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) +
      scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) +
      labs(y = "Gene", x= "Neuron") + theme(panel.grid = element_line(size = 0.5, color = "grey85"))
    
  fnh <-
    function() {
      withProgress(message = "Generating heatmap Plot...", value = 0, { g })
    }
  }

  output$downloadheatmap <-
    downloadHandler(
      filename = function() {
        filename = paste("Heatmap.png", sep = "")
      },
      content = function(file) {
        ggsave(
          fnh(),
          file = file,
          height = 200,
          width = 500,
          units = "mm" ,
          limitsize = FALSE,
          device = "png"
        )
      }
    )

  output$heatmap <- renderPlot(g)
  
  output$dynamic <- renderUI({
    #req(input$plot_hover)
    verbatimTextOutput("vals", placeholder = TRUE)
  })
  
  output$vals <- renderPrint({
    hover <- input$plot_hover 
    #print(input$plot_hover) # list
    y <- nearPoints(flp.neuron.scaled, input$plot_hover)
    req(nrow(y) != 0)
    y
  })
})

### From file

  observeEvent(input$PlotHeatmap2, {
    ds <- input$dataset_heatmap
    inFile <- input$file1
    ss<-read.table(inFile$datapath, header=FALSE)$V1
    
    ss <- strsplit(as.character(input$genelist), "\n| |\\,|\t")
    ss <- as.data.frame(ss)[,1]
    ss <- unique(c(ss, filter(gene_list, gene_id %in% ss | seqnames %in% ss)$gene_name))
    
    mis <- ss[ss %in% c(gene_list$gene_id, gene_list$gene_name, gene_list$seqnames) & !ss %in% med.scaled.long$gene_name & ss %in% gene_list$gene_name]
    mis_all <- ss[ss %in% c(gene_list$gene_id, gene_list$gene_name, gene_list$seqnames) & !ss %in% L4.TPM.raw.scaled.long$gene_name & ss %in% gene_list$gene_name]
    
    if(ds=="Neurons only"){
      L4.TPM=L4.TPM.medium
      heatmapdata=med.scaled.long
      cc = colnames(ths)[-c(1,130,131)]
      missing = mis} else {
        L4.TPM=as(L4.all.TPM.raw,"dgCMatrix")
        heatmapdata=L4.TPM.raw.scaled.long
        cc=colnames(L4.all.TPM.raw)
        missing = mis_all
      }
    
    head(heatmapdata)
    print(ds)
    
    flp.neuron.scaled <- heatmapdata[which(heatmapdata$gene_name %in% ss),]
    flp.ids <- as.character(vlookup(unique(flp.neuron.scaled$gene_name), gene_list, result_column = 1, lookup_column = 2))
    flp.expr <- L4.TPM[flp.ids,, drop=FALSE]
    if ( nrow(flp.expr) >1 ) {
      flp.neuron.order <- pheatmap(flp.expr, scale = "row")
      flp.neuron.order <- flp.neuron.order[["tree_row"]]$order
      flp.neuron.order <- flp.ids[flp.neuron.order]
    } else {
      flp.neuron.order <- flp.ids[1]
    }
    flp.neuron.order <- as.character(vlookup(flp.neuron.order, gene_list))
    flp.neuron.scaled$gene_name <- factor(flp.neuron.scaled$gene_name, levels = c(rev(flp.neuron.order), missing))
    #flp.neuron.scaled$gene_name <- fct_rev(flp.neuron.scaled$gene_name)
    
    for( i in missing ){
      dff <- data.frame(gene_name=i, cell.type= cc, scaled.expr=0, prop=0, Modality="NA")
      if(ds!="Neurons only"){colnames(dff)[5]<-"tissue"}
      flp.neuron.scaled <- rbind(flp.neuron.scaled, dff)
    }
    
    if (ds!="Neurons only"){
      
      if(nrow(flp.expr) >1){
        g <- ggplot(flp.neuron.scaled, aes(y = gene_name, x = cell.type)) + 
          geom_point(aes(color = scaled.expr, size = prop)) +
          theme(panel.background = element_blank(), axis.text.y.left = element_text(size = 12, color = "black"),
                legend.key = element_blank(),
                legend.text = element_text(color = "black", size = 14),
                legend.title = element_text(color = "black", size = 16),
                axis.title = element_text(size = 16, color = "black")) + 
          scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) + 
          scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) + 
          labs(y = "Gene", x = "Tissue") + theme(panel.grid = element_line(size = 0.2, color = "grey85")) +
          facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") + 
          theme(strip.placement = "outside", 
                strip.background.x = element_blank(),
                axis.text.x.bottom = element_blank(),
                strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) 
      } else {
        
        g <- ggplot(flp.neuron.scaled, aes(y = gene_name, x = cell.type)) + 
          geom_point(aes(color = scaled.expr, size = prop)) +
          theme(panel.background = element_blank(), axis.text.y.left = element_text(size = 12, color = "black"),
                legend.key = element_blank(),
                legend.text = element_text(color = "black", size = 14),
                legend.title = element_text(color = "black", size = 16),
                axis.title = element_text(size = 16, color = "black")) + 
          scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) + 
          scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) + 
          labs(y = "Gene", x = "Tissue") + theme(panel.grid = element_line(size = 0.2, color = "grey85"))
      }
      
      pg <- ggplotGrob(g)
      
      for(i in which(grepl("strip-b", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
      }
      
      fnh <-
        function() {
          withProgress(message = "Generating heatmap Plot...", value = 0, {
            grid::grid.draw(pg)      })
        }
    } else {
      
      g<-ggplot(flp.neuron.scaled, aes(y =gene_name, x = cell.type)) +
        geom_point(aes(color = scaled.expr, size = prop)) +
        #geom_point(data = NULL, aes(y = vec), pch = NA) +
        theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 7, color = "black", hjust = 1),
              panel.background = element_blank(), axis.text.y.left = element_text(size = 7, color = "black")) +
        scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) +
        scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) +
        labs(y = "Gene", x= "Neuron") + theme(panel.grid = element_line(size = 0.5, color = "grey85"))
      
      fnh <-
        function() {
          withProgress(message = "Generating heatmap Plot...", value = 0, { g })
        }
    }
    
    output$downloadheatmap <-
      downloadHandler(
        filename = function() {
          filename = paste("Heatmap.png", sep = "")
        },
        content = function(file) {
          ggsave(
            fnh(),
            file = file,
            height = 200,
            width = 500,
            units = "mm" ,
            limitsize = FALSE,
            device = "png"
          )
        }
      )
    
    output$heatmap <- renderPlot(g)
    
    output$dynamic <- renderUI({
      #req(input$plot_hover)
      verbatimTextOutput("vals", placeholder = TRUE)
    })
    
    output$vals <- renderPrint({
      hover <- input$plot_hover 
      #print(input$plot_hover) # list
      y <- nearPoints(flp.neuron.scaled, input$plot_hover)
      req(nrow(y) != 0)
      y
    })
  })

}