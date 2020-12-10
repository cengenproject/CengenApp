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
library(shinybusy)
library(tidyverse)
library(pheatmap)
library(expss)
library(plotly)
#library(cairo)
#library(MAST)
options(repos = BiocManager::repositories())
options("repos")
options(shiny.maxRequestSize = 400 * 1024 ^ 2)
source("Functions.R")

options(scipen = 0)
options(digits = 2)

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
    "WBGene00001555"
  )



## UI ----
ui <- fluidPage(
  theme = "Theme.min.css",
  #shinytheme("flatly"),
  #background-color: #2f2d2d;
  #shinythemes::themeSelector(),
  tags$head(tags$style(
    HTML(".shiny-output-error-validation {color: red;}")
  )),
  
  # App title ----
  titlePanel(
    "Discovery and analysis of the C. elegans Neuronal Gene Expression Network -- CeNGEN"
  ),
  hr(),
  add_busy_spinner(
    spin = "double-bounce",
    color = "orange",
    timeout = 100,
    position = c("top-right"),
    onstart = TRUE,
    margins = c(10, 10),
    height = "50px",
    width = "50px"
  ),
  #add_busy_bar(color = "darkorange"),
  # Side panel entering what to plot ----
  
  
  # Main panel showing plots ----
  tabsetPanel(
    type = "pills",
    id = "tabs",
    ### Cell types Panel ----
    tabPanel(
      "Gene expression by cell type",
      fluidPage(
        hr(),
        h6(
          "Find which genes are expressed in each cell type ordered by expression level or which cell type express a particular gene."
        ),
        h6("Select the threshold of expression or unfiltered data"),
        h6("4 More stringent, 1 less stringent"),
        hr(),
        fluidRow(
          column(1),
          column(
            4,
            selectInput(
              inputId = "Tcell_name",
              label = "Select cell type",
              choices = colnames(ths)[2:129],
              selected = "ADA"
            ),
            selectInput(
              inputId = "Tcell_cut",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered"),
              selected = 2
            ),
            actionButton("TCell", "Expressed genes", icon = icon("hand-o-right"))
            
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            textInput(
              inputId = "Tgene_name",
              label = "Type gene name",
              value = "zig-4"
            ),
            selectInput(
              inputId = "Tgene_cut",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered"),
              selected = 2
            ),
            actionButton("TGene", "Which cell types", icon = icon("hand-o-right"))
          ),
          column(
            2,
            offset = 0,
            style = 'padding:5px;',
            textInput(
              inputId = "Tgene_name_batch",
              label = "Query multiple genes for download",
              value = "zig-4,ins-26"
            ),
            selectInput(
              inputId = "Tgene_cut_batch",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered"),
              selected = 2
            ),
            downloadButton("TGeneBatch", "Download batch"),
            span(textOutput("textb"), style =
                   "color:red")
            
          )
        ),
        br(),
        fluidRow(
          column(1),
          column(
            4,
            br(),
            DT::dataTableOutput("Tcell_name_table"),
            br(),
            uiOutput("get_download_gene")
            #downloadButton('downloadGene', "Download table")
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            br(),
            DT::dataTableOutput("Tgene_name_table"),
            br(),
            #downloadButton('downloadCell', "Download table"),
            uiOutput("get_download_cell"),
            span(textOutput("text1"), style =
                   "color:red")
          )
        )
      )
    ),
    
    tabPanel(
      "Find markers based on percentage of expression",
      fluidPage(
        hr(),
        h6(
          "Find which genes are expressed in a group of cell types and not in other cell types."
        ),
        fluidRow(
          #column(1),
          column(
            4,
            textInput(
              inputId = "String1",
              label = "Group 1",
              value = "AWC_ON,AWC_OFF"
            ),
            textInput(
              inputId = "Expressed",
              label = "Minimum percentage of cells expressing the gene",
              value = 65
            ),
            actionButton("Filter", "Run query", icon = icon("hand-o-right"))
            
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            textInput(
              inputId = "String2",
              label = "Group 2",
              value = "SMD,SIB"
            ),
            textInput(
              inputId = "NonExpressed",
              label = "Maximum percentage of cells expressing the gene",
              value = 2
            ),
            downloadButton('downloadQuery', "Download table")
          )
        ),
        fluidRow(
          #column(1),
          column(4,
                 br(), DT::dataTableOutput("YesExpressed")),
          column(4,
                 br(),
                 DT::dataTableOutput("NoExpressed")),
          
          column(4,
                 br(),
                 DT::dataTableOutput("Result")
                 #span(textOutput("text3"), style="color:red"))
          )
          
        )
      )
    ),
      ### Enriched types Panel ----
      tabPanel(
        "Enriched Genes by Cell Type",
        fluidPage(
          hr(),
          h6(
            "Find genes overexpressed in one cell type compared to the rest of cells in the dataset (neurons only or all cell types)."
          ),
          textOutput("TopMarkers cell plot"),
          selectInput(
            inputId = "dataset2",
            label = "Choose dataset",
            choices = c("All cell types", "Neurons only")
          ),
          conditionalPanel(
            "input.dataset2 == 'Neurons only'",
            selectInput(
              inputId = "Markers",
              label = "Select cluster",
              choices = sort(unique(allNeurons$`Neuron.type`)),
              selected = "ADA"
            ),
            textInput(
              inputId = "top",
              label = "Show top X genes",
              value = "30"
            ),
            DT::dataTableOutput("MarkTable"),
            downloadButton('downloadMarkers', "Download table"),
            h6("HEADER LEGEND:"),
            h6(
              "p-val and p_val_adj: nominal and adjusted P-values of the test, respectively."
            ),
            h6(
              "pct.1 and pct.2: The percentage of cells where the gene is detected in the first or second group"
            ),
            h6(
              "avg_logFC: Log of the expression fold change between group 1 and group 2."
            )
          ),
          conditionalPanel(
            "input.dataset2 == 'All cell types'",
            selectInput(
              inputId = "Markers2",
              label = "Select cluster",
              choices = sort(unique(allCells$`Neuron.type`)),
              selected = "Sperm"
            ),
            textInput(
              inputId = "top2",
              label = "Show top X genes",
              value = "30"
            ),
            DT::dataTableOutput("MarkTable2"),
            downloadButton('downloadMarkers2', "Download table"),
            h6("HEADER LEGEND:"),
            h6(
              "p-val and p_val_adj: nominal and adjusted P-values of the test, respectively."
            ),
            h6(
              "pct.1 and pct.2: The percentage of cells where the gene is detected in the first or second group"
            ),
            h6(
              "avg_logFC: Log of the expression fold change between group 1 and group 2."
            )
          )
        )
      ),
      ### DEX panel ----
      tabPanel(
        "Find Differential Expression between Cell Types",
        fluidPage(
          hr(),
          h6("The calculation can take a few minutes."),
          hr(),
          fluidRow(
            column(
              3,
              selectInput(
                inputId = "ClusterCells1",
                label = "Select cluster 1",
                choices = sort(unique(allCells$`Neuron.type`)),
                selected = "Sperm"
              ),
              selectInput(
                inputId = "ClusterCells2",
                label = "Select cluster 2",
                choices = sort(unique(allCells$`Neuron.type`)),
                selected = "Intestine"
              ),
              actionButton("DEXButton", "Calculate DEX", icon = icon("hand-o-right"))
              
            ),
            column(
              3,
              selectInput(
                inputId = "test",
                label = "Select statistical test",
                choices = c("wilcox", "bimod", "roc", "t", "LR")
              ),
              textInput(
                inputId = "topM2",
                label = "Show top X genes",
                value = "30"
              )
              
            ),
            column(
              3,
              textInput(
                inputId = "batch1",
                label = "Group 1: Introduce cell types separated with ,",
                value = "AVM,AWA,AWB"
              ),
              textInput(
                inputId = "batch2",
                label = "Group 2: Introduce cell types, NEURONS or ALL",
                value = "AWC_ON,AWC_OFF"
              ),
              actionButton(
                "DEXButtonBatch",
                "Calculate DEX in groups",
                icon = icon("hand-o-right")
              )
            )
          ),
          br(),
          h6("HEADER LEGEND:"),
          h6(
            "p-val and p_val_adj: nominal and adjusted P-values of the test, respectively."
          ),
          h6(
            "pct.1 and pct.2: The percentage of cells where the gene is detected in the first or second group"
          ),
          h6(
            "avg_logFC: Log of the expression fold change between group 1 and group 2."
          ),
          br(),
          
          DT::dataTableOutput("MarkTable_ClusterCells"),
          span(textOutput("text2"), style =
                 "color:red"),
          br(),
          DT::dataTableOutput("MarkTable_Batch"),
          downloadButton('downloadDEX', "Download table")
        )
      ),
      ### Single cell panel ----
      tabPanel("Single cell plot", fluidPage(
        hr(),
        column(
          3,
          h4('Integrated single-cell (10X) analysis'),
          wellPanel(
            selectInput(
              inputId = "dataset",
              label = "Choose dataset",
              choices = c("All cell types", "Neurons only")
            ),
            selectInput(
              inputId = "Plots",
              label = "Features to plot",
              choices = c("Metadata", "Individual genes", "Colocalization"),
              selected = "Individual genes"
            ),
            # One feature ----
            conditionalPanel(
              ### METADATA PANELS neurons ----
              condition = "input.Plots == 'Metadata' && input.dataset == 'Neurons only'",
              selectInput(
                "featuresMetadata",
                "Metadata",
                colnames(allNeurons@meta.data)
              ),
              
              ### CLUSTERS PANELS
              conditionalPanel(
                condition = "input.featuresMetadata == 'Neuron.type' && input.dataset == 'Neurons only'",
                selectInput(
                  "Cluster",
                  "Highlight cells",
                  choices = unique(allNeurons$Neuron.type)
                )
              ),
              conditionalPanel(
                condition = "input.featuresMetadata == 'Experiment' && input.dataset == 'Neurons only'",
                selectInput(
                  "Cluster2",
                  "Highlight cells",
                  choices = unique(allNeurons$Experiment)
                )
              ),
              conditionalPanel(
                condition = "input.featuresMetadata == 'Detection' && input.dataset == 'Neurons only'",
                selectInput("Cluster3", "Highlight cells", choices = unique(allNeurons$Detection))
              )
            ),
            
            conditionalPanel(
              ### METADATA PANELS all cells ----
              condition = "input.Plots == 'Metadata' && input.dataset == 'All cell types'",
              selectInput(
                "featuresMetadata2",
                "Metadata",
                colnames(allCells@meta.data)
              ),
              ### CLUSTERS PANELS
              conditionalPanel(
                condition = "input.featuresMetadata2 == 'Neuron.type' && input.dataset == 'All cell types'",
                selectInput(
                  "Cluster.2",
                  "Highlight cells",
                  choices = unique(allCells$Neuron.type)
                )
              ),
              conditionalPanel(
                condition = "input.featuresMetadata2 == 'Experiment' && input.dataset == 'All cell types'",
                selectInput(
                  "Cluster2.2",
                  "Highlight cells",
                  choices = unique(allCells$Experiment)
                )
              ),
              conditionalPanel(
                condition = "input.featuresMetadata2 == 'Detection' && input.dataset == 'All cell types'",
                selectInput("Cluster3.2", "Highlight cells", choices = unique(allCells$Detection))
              ),
              conditionalPanel(
                condition = "input.featuresMetadata2 == 'Tissue.type' && input.dataset == 'All cell types'",
                selectInput("Cluster4", "Highlight cells", choices = unique(allCells$Tissue.type))
              )
            ),
            # One gene ----
            conditionalPanel(
              condition = "input.Plots == 'Individual genes'",
              textInput(
                inputId = "GeneName",
                label = "Gene name (Symbol)",
                value = "ric-4"
              )
            ),
            # Two gene colocalization ----
            conditionalPanel(
              condition = "input.Plots == 'Colocalization'",
              textInput(
                inputId = "ColocalizationGenes1",
                label = "Gene 1",
                value = "zig-4"
              ),
              hr(),
              textInput(
                inputId = "ColocalizationGenes2",
                label = "Gene 2",
                value = "F23D12.3"
              )
            ),
            fluidRow(column(
              6, actionButton("PlotButton", "Plot data", icon = icon("hand-o-right"))
            ))
          )
        ),
        
        column(
          9,
          textOutput("Single cell plot"),
          conditionalPanel(
            condition = "input.Plots == 'Metadata'",
            downloadLink("downloadClust", "Download Plot"),
            plotOutput(
              "SinglePlot",
              width = "100%",
              dblclick = "plotPermanent_dblclick",
              brush = brushOpts(id = "plotPermanent_brush", resetOnNew = TRUE)
            ),
            h5(
              "Select region and do double-click to zoom in and double-click to zoom out."
            ),
            hr(),
            downloadLink("downloadHigh", "Download Plot"),
            plotOutput("High", width = "100%")
          ),
          conditionalPanel(
            condition = "input.Plots == 'Colocalization'",
            sliderInput(
              "blend",
              "Choose blend",
              min = 0,
              max = 1,
              value = 0.25
            ),
            plotOutput("SinglePlotDouble", width =
                         "50%"),
            downloadLink("downloadCol", "Download plot"),
            h4(span(textOutput("gene1utr"), style =
                      "color:red")),
            h4(span(textOutput("gene2utr"), style =
                      "color:red")),
            hr(),
            h5(
              "Select region in the panel below and do double-click to zoom in and double-click to zoom out."
            ),
            plotOutput(
              "SinglePlotPermanent",
              width = "50%",
              dblclick = "plotPermanent_dblclick",
              brush = brushOpts(id = "plotPermanent_brush", resetOnNew = TRUE)
            ),
            hr()
          ),
          conditionalPanel(
            condition = "input.Plots == 'Individual genes'",
            h5(
              "Select region and do double-click to zoom in and double-click to zoom out."
            ),
            plotOutput(
              "SinglePlot2",
              width = "100%",
              dblclick = "plotPermanent_dblclick",
              brush = brushOpts(id = "plotPermanent_brush", resetOnNew = TRUE)
            ),
            h4(span(textOutput("geneUTR"), style =
                      "color:red")),
            hr(),
            downloadLink("downloadExp", "Download plot"),
            plotOutput("FeaturePlot", width = "100%"),
            hr(),
            downloadLink("downloadViolin", "Download plot"),
            plotOutput("Violin", width = "100%")
          )
        )
      )),
      tabPanel(
        "Heatmaps of gene expression",
        fluidPage(
          hr(),
          h6(
            "Introduce a list of genes to display a heatmap showing expression levels and proportion of cells expressing the genes."
          ),
          textAreaInput(
            inputId = "genelist",
            label = "Introduce a list of genes",
            value = "flp-1\nflp-2,flp-3,WBGene00001447\nWBGene00001448\nflp-6\nflp-7\nflp-8\nflp-9\nflp-10\nflp-11\nflp-12\nflp-13\nflp-14\nflp-15\nflp-16\nflp-17\nflp-18\nflp-19\nflp-20\nflp-21\nflp-22\nflp-23\nflp-24\nflp-25\nflp-26\nflp-27\nflp-28\nflp-32\nflp-33\nflp-34",
            width = "500px",
            height = "100px"
          ),
          actionButton(
            "PlotHeatmap",
            "Plot heatmap from list",
            icon = icon("hand-o-right")
          ),
          hr(),
          fluidRow(
            column(3,fileInput("file1", NULL,
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      "txt")
            )),
            #column(3,actionButton("resetFile", "Clear uploaded file")),
            column(3,   actionButton(
              "PlotHeatmap2",
              "Plot heatmap from file",
              icon = icon("hand-o-right")
            ))
          ),
        
          hr(),
          h6(
            "You can identify circles by clicking on them."
          ),
         
          
          uiOutput("dynamic"),
          hr(),
          plotOutput("heatmap", width = "100%", hover = "plot_hover"),
          
          hr(),
          downloadLink("downloadheatmap", "Download plot")
          
            
        )
      )
    )
  )
 
