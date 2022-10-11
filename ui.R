#### CenGenAPP 2018-2020
#### Report bugs to: gabrielsantperebaro@gmail.com

require(shiny)
library(shinyjs)
library(shinythemes)

library(DT)
library(tidyverse)
library(pheatmap)

library(ggridges)
library(cowplot)
library(shinybusy)
library(plotly)

library(expss)
library(xlsx)

library(edgeR)
library(qs)
library(Matrix)
library(limma)

options(repos = BiocManager::repositories())
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
    "WBGene00001555",
    "WBGene00206533",
    "WBGene00011964",
    "WBGene00018172",
    "WBGene00016259",
    "WBGene00023407"
  )

msg <- filter(gene_list, gene_id %in% utr)$gene_name %>% paste(., collapse = ", ")


## UI ----
ui <- fluidPage(
  tags$head(includeHTML(("google-analytics-script2.html"))),
  theme = "Theme.min.css",
  tags$head(tags$style(
    HTML(".shiny-output-error-validation {color: red;}")
  )),
  
  # App title ----
  titlePanel(
    "Discovery and analysis of the C. elegans Neuronal Gene Expression Network -- CeNGEN"
  ),
  tags$a(href="http://www.cengen.org/single-cell-rna-seq/", "CeNGENApp Help and Documentation"),
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
  
  # Main panel showing plots ----
  tabsetPanel(
    type = "pills",
    id = "tabs",
    ### Gene expression by cell type Panel ----
    tabPanel(
      "Gene expression by cell type",
      fluidPage(
        hr(),
        h6(
          "Find all genes expressed in a given cell type, or all cell types expressing a given gene (or group of genes)."
        ),
        h6("Select one of four thresholds for expression:"),
        h6("1 (least stringent) to 4 (most stringent) or select unfiltered data"),
        h6("Choose All Cells Unfiltered to query the entire unfiltered dataset, including non-neuronal cells"),
        hr(),
        fluidRow(
          column(1),
          column(
            4,
            selectInput(
              inputId = "Tcell_name",
              label = "Select cell type",
              choices = all_cell_types,
              selected = "ADA"
            ),
            selectInput(
              inputId = "Tcell_cut",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered", "All Cells Unfiltered"),
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
              choices = c(1:4, "Unfiltered", "All Cells Unfiltered"),
              selected = 2
            ),
            actionButton("TGene", "Which cell types", icon = icon("hand-o-right"))
          ),
          column(
            2,
            offset = 0,
            style = 'padding:5px;',
            textAreaInput(
              inputId = "Tgene_name_batch",
              label = "Query multiple genes for download",
              value = "flp*\nWBGene00001447\nWBGene00001448,zig-4",
              width = "500px",
              height = "100px"
            ),
            selectInput(
              inputId = "Tgene_cut_batch",
              label = "Select threshold",
              choices = c(1:4, "Unfiltered", "All Cells Unfiltered" ),
              selected = 2
            ),
            downloadButton("TGeneBatch", "Download batch"),
            span(textOutput("textb"), style =
                   "color:red")
            
          )
        ),
        br(),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        fluidRow(
          column(1),
          column(
            4,
            br(),
            span(textOutput("Error1"), style ="color:red"),
            DT::dataTableOutput("Tcell_name_table"),
            br(),
            uiOutput("get_download_gene")
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
    
    ### Find markers based on percentage of expression Panel ----
    tabPanel(
      "Find markers based on percentage of expression",
      fluidPage(
        hr(),
        h6(
          "Find genes expressed in one group of cell types and not in another group based on the percentages of cells expressing the gene."
        ),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        br(),
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
    ### Enriched Genes by cell type Panel ----
    tabPanel(
      "Enriched Genes by cell type",
      fluidPage(
        hr(),
        h6(
          "Find genes differentially expressed in one cell type compared to all other cells in the dataset (neurons only or all cell types). 
            This is NOT a comprehensive list of genes detected in each cell type."
        ),
        h6("Please see Gene Expression by Cell type for a comprehensive list of expression values of each gene in a given cell type."),
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
            choices = all_neuron_types,
            selected = "SIA"
          ),
          textInput(
            inputId = "top",
            label = "Show top X genes",
            value = "100"
          ),
          h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
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
            choices = all_cell_types,
            selected = "SIA"
          ),
          textInput(
            inputId = "top2",
            label = "Show top X genes",
            value = "100"
          ),
          h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
          
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
    
    
    ### Find Differential Expression between Cell Types Panel ----
    tabPanel(
      "Find Differential Expression between Cell Types",
      fluidPage(
        hr(),
        h6("Find differentially expressed genes between two cell types or two groups of cell types."),
        h6("The calculation can take a few minutes."),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        hr(),
        fluidRow(
          column(
            4,
            selectInput(
              inputId = "batch1",
              label = "Select Group 1",
              choices = all_cell_types,
              selected = "AVL",
              multiple = TRUE
            ),
            selectInput(
              inputId = "batch2",
              label = "Select Group 2: Introduce cell types, NEURONS or ALL",
              choices = c("ALL","NEURONS", all_cell_types),
              selected = c("RME_DV","RME_LR"),
              multiple = TRUE
            ),
            actionButton("DEXButton", "Calculate DEX", icon = icon("hand-o-right"))
            
          ),
          column(
            3,
            selectInput(
              inputId = "test",
              label = "Select statistical test",
              choices = c("Wilcoxon on single cells", "Pseudobulk: Wilcoxon", "Pseudobulk: edgeR pairwise exact test")
            ),
            textInput(
              inputId = "topM2",
              label = "Show top X genes",
              value = "100"
            )
            
          )
        ),
        br(),
        htmlOutput("pseudobulk_metadata", container = h6),
        br(),
        
        span(textOutput("text_error_dex"), style =
               "color:red"),
        br(),
        DT::dataTableOutput("MarkTable_Batch"),
        downloadButton('downloadDEX', "Download table"),
        htmlOutput("legend_de_columns", container = h6)
      )
    ),
    
    
    ### Heatmaps of gene expression Panel ----
    tabPanel(
      "Heatmaps of gene expression",
      fluidPage(
        hr(),
        h6(
          "Display a heatmap showing relative expression and proportion of cells expressing a gene or group of genes across all neurons. This function uses data from threshold 2. Color shows relative scaled expression for each gene across neuron types, and is not comparable between genes."
        ),
        h6( paste0("WARNING: Expression values for ",msg," are unreliable as they have been overexpressed to generate transgenic strains."), style="color:orange"),
        textAreaInput(
          inputId = "genelist",
          label = "Introduce a list of genes",
          value = "flp*\nWBGene00001447\nWBGene00001448,zig-4",
          width = "500px",
          height = "100px"
        ),
        
        selectInput(
          inputId = "dataset_heatmap",
          label = "Choose dataset: Neurons (threshold 2), All cells (unfiltered)",
          choices = c("Neurons only", "All cell types")
        ),    
        
        
        actionButton(
          "PlotHeatmap","Plot heatmap from list",
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
        
        h6(
          "You can identify circles by clicking on them."
        ),
        
        div(style="height:30px;width:800px;padding-left:10px;padding-right:10px;background-color:#ffffff;",fluidRow(verbatimTextOutput("vals", placeholder = TRUE))),
        #uiOutput("dynamic"),
        br(),
        br(),
        br(),
        plotOutput("heatmap", width = "100%",hover = "plot_hover"),
        
        hr(),
        downloadLink("downloadheatmap", "Download plot")
        
        
      )
    )
  )
)

