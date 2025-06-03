#### CeNGEN-app 2018-2022
#### Report bugs at https://github.com/cengenproject/CengenApp/issues


cat("############# Loading libraries #################\n")

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
library(svglite)

library(expss)

library(BiocManager)
library(edgeR)
library(qs)
library(Matrix)
library(limma)

cat("############## Setting options ##################\n")

options(repos = BiocManager::repositories())
options(shiny.maxRequestSize = 400 * 1024 ^ 2)

options(scipen = 0)
options(digits = 2)

cat("############ End initializations ################\n")


source("Functions.R")



warning_msg <- if(!all(is.na(unreliable_gene_ids))){
  p( paste0(
    "WARNING: Expression values for ",
    filter(gene_list, gene_id %in% unreliable_gene_ids)$gene_name |> paste(collapse = ", "),
    " are unreliable as they have been overexpressed to generate transgenic strains."
  ),
  style="color:orange")
} else{
  ""
}


## UI ----
ui <- fluidPage(
  
  #~ HTML head ----
  
  
  tags$head(
    
    includeHTML(("www/google-analytics-script2.html")),
    
    tags$style(
      HTML(
        ".shiny-output-error-validation {color: red;}",
        "html {
          font-size: 14px !important;
        }",
        "body {
          padding-bottom: 100px;
        }"
      )
    ),
    
    tags$link(rel="shortcut icon", href=favicon),
    
  ),
  
  theme = "startbs_landing_styles.css",
  
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
  
  #~ App header ----
  div(
    tags$img(src = icon_big, style = "margin: 25px"),
    div(
      titlePanel(
        paste(dataset, "CeNGEN")
      ),
      p("Discovery and analysis of the C. elegans Neuronal Gene Expression Network", 
        style = "font-size: 1.2em; color: #555; margin-top: -10px; margin-bottom: 0;"),
      style = "flex: 1;"
    ),
    style = paste0(
      "display: flex; align-items: center;  margin-bottom: 30px; padding: 20px;",
      "background: linear-gradient(to bottom, ",
      bg_color, " 0%, ",
      bg_color, " 40%, ",
      "white 100%);"
    )
  ),
  
  div(
    "This app enables analysis of the ",
    
    tags$img(src = paste0(dataset, ".png")),
    dataset, " dataset. Other datasets: ",
    HTML(paste0('<img src="',other_apps$name,'.png"/> <a href="https://cengen.shinyapps.io/',other_apps$url,'">',other_apps$name,'</a>',
                collapse = ", ")),
    tags$p("Additional details and documentation ",
           tags$a("on the cengen.org page",
                  href = "https://www.cengen.org/single-cell-rna-seq/"),
           "."),
    class = "alert alert-primary"
    
  ),
  hr(),
  
  # Main panel showing plots ----
  tabsetPanel(
    type = "pills",
    id = "tabs",
    ### Gene expression by cell type Panel ----
    tabPanel(
      "Gene expression by cell type",
      fluidPage(
        hr(),
        p(
          "Find all genes expressed in a given cell type, or all cell types expressing a given gene (or group of genes).",
          br(),
          "Select one of four thresholds for expression:\n",
          br(),
          "1 (least stringent) to 4 (most stringent) or select unfiltered data\n",
          br(),
          "Choose All Cells Unfiltered to query the entire unfiltered dataset, including non-neuronal cells"
        ),
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
            actionButton("TCell", "Expressed genes", class = "btn-info")
            
          ),
          
          
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
            actionButton("TGene", "Which cell types", class = "btn-info")
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
            downloadButton("TGeneBatch", "Download batch",
                           class = "btn-secondary"),
            span(textOutput("textb"), style =
                   "color:red")
            
          )
        ),
        br(),
        warning_msg,
        fluidRow(
          column(1),
          column(
            4,
            br(),
            span(textOutput("Error1"), style ="color:red"),
            DT::dataTableOutput("Tcell_name_table"),
            br(),
            uiOutput("get_download_gene",
                     class = "btn-secondary")
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            br(),
            DT::dataTableOutput("Tgene_name_table"),
            br(),
            #downloadButton('downloadCell', "Download table"),
            uiOutput("get_download_cell",
                     class = "btn-secondary"),
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
        p(
          "Find genes expressed in one group of cell types and not in another group based on the percentages of cells expressing the gene."
        ),
        warning_msg,
        br(),
        fluidRow(
          #column(1),
          column(
            4,
            textInput(
              inputId = "PCT_group1",
              label = "Group 1",
              value = "AWC_ON,AWC_OFF"
            ),
            textInput(
              inputId = "PCT_expressed_threshold",
              label = "Minimum percentage of cells expressing the gene",
              value = 65
            ),
            actionButton("Filter", "Run query", class = "btn-info")
            
          ),
          #column(width = 1, offset = 0, style='padding:5px;'),
          column(
            4,
            textInput(
              inputId = "PCT_group2",
              label = "Group 2",
              value = "SMD,SIB"
            ),
            textInput(
              inputId = "PCT_notExpressed_threshold",
              label = "Maximum percentage of cells expressing the gene",
              value = 2
            ),
            downloadButton('downloadQuery', "Download table",
                           class = "btn-secondary")
          )
        ),
        fluidRow(
          #column(1),
          column(4,
                 div(
                   h5("Expressed in group 1", class = "card-title"),
                   DT::dataTableOutput("YesExpressed"),
                   class = "card-body",
                   style = "overflow-x: auto;"
                 ),
                 class = "card"),
          column(4,
                 div(
                   h5("Not expressed in group 2", class = "card-title"),
                   DT::dataTableOutput("NoExpressed"),
                   class = "card-body",
                   style = "overflow-x: auto;"
                 ),
                 class = "card"),
          
          column(4,
                 div(
                   h5("Result (expressed in group 1, not group 2)", class = "card-title"),
                   DT::dataTableOutput("Result"),
                   class = "card-body"
                 ),
                 class = "card"),
          class = "mt-3"
        )
      )
    ),
    ### Enriched Genes by cell type Panel ----
    tabPanel(
      "Enriched Genes by cell type",
      fluidPage(
        hr(),
        p(
          "Find genes differentially expressed in one cell type compared to all other cells in the dataset (neurons only or all cell types).",
          "This is NOT a comprehensive list of genes detected in each cell type.",
          br(),
          "Please see Gene Expression by Cell type for a comprehensive list of expression values of each gene in a given cell type."
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
            choices = all_neuron_types,
            selected = "SIA"
          ),
          textInput(
            inputId = "top",
            label = "Show top X genes",
            value = "100"
          ),
          warning_msg,
          DT::dataTableOutput("MarkTable"),
          downloadButton('downloadMarkers', "Download table",
                         class = "btn-secondary"),
          p("HEADER LEGEND:",
            br(),
            "p-val and p_val_adj: nominal and adjusted P-values of the test, respectively.",
            br(),
            "pct.1 and pct.2: The percentage of cells where the gene is detected in the first or second group",
            br(),
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
          warning_msg,
          
          DT::dataTableOutput("MarkTable2"),
          downloadButton('downloadMarkers2', "Download table",
                         class = "btn-secondary"),
          p("HEADER LEGEND:",
            br(),
            "p-val and p_val_adj: nominal and adjusted P-values of the test, respectively.",
            br(),
            "pct.1 and pct.2: The percentage of cells where the gene is detected in the first or second group",
            br(),
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
        p("Find differentially expressed genes between two cell types or two groups of cell types.",
          br(),
          "Note, this computation is performed on demand. Comparisons of large number of cells can take several minutes and lead to app disconnections. If this becomes a problem, consider using a Pseudobulk test or running a local version of the app."
        ),
        warning_msg,
        hr(),
        fluidRow(
          column(
            4,
            selectInput(
              inputId = "DEgroup1",
              label = "Select Group 1",
              choices = all_cell_types,
              selected = "AVL",
              multiple = TRUE
            ),
            selectInput(
              inputId = "DEgroup2",
              label = "Select Group 2: Introduce cell types, NEURONS or ALL",
              choices = c("ALL","NEURONS", all_cell_types),
              selected = c("RME_DV","RME_LR"),
              multiple = TRUE
            ),
            actionButton("DEbutton", "Calculate DEX", class = "btn-info")
            
          ),
          column(
            3,
            selectInput(
              inputId = "DEtest",
              label = "Select statistical test",
              choices = c("Wilcoxon on single cells", "Pseudobulk: Wilcoxon", "Pseudobulk: edgeR pairwise exact test")
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
        downloadButton('downloadDEX', "Download table",
                       class = "btn-secondary"),
        htmlOutput("legend_de_columns", container = h6)
      )
    ),
    
    
    ###Sex DE ----
    
    if(dataset == "male"){
      
      tabPanel(
        "Male vs hermaphrodite DE",
        fluidPage(
          hr(),
          p(
            "Find differentially expressed genes between male and hermaphrodite cell types or two groups of cell types.",
            br(),
            "Note, this computation is performed on demand. Comparisons of large number of cells can take several minutes and lead to app disconnections. If this becomes a problem, consider using a Pseudobulk test or running a local version of the app."
          ),
          warning_msg,
          hr(),
          fluidRow(
            column(
              4,
              selectInput(
                inputId = "SDEmale",
                label = "Select male cells",
                choices = male_all_cell_types,
                selected = "ADF",
                multiple = TRUE
              ),
              selectInput(
                inputId = "SDEherm",
                label = "Select hermaphrodite cells",
                choices = herm_all_cell_types,
                selected = c("ADF"),
                multiple = TRUE
              ),
              actionButton("SDEbutton", "Calculate DEX", class = "btn-info")
              
            ),
            column(
              3,
              selectInput(
                inputId = "SDEtest",
                label = "Select statistical test",
                choices = c("Wilcoxon on single cells", "Pseudobulk: Wilcoxon", "Pseudobulk: edgeR pairwise exact test"),
                selected = "Pseudobulk: edgeR pairwise exact test"
              )
              
            )
          ),
          br(),
          htmlOutput("SDEpseudobulk_metadata", container = h6),
          br(),
          
          span(textOutput("SDEtext_error_dex"), style =
                 "color:red"),
          br(),
          DT::dataTableOutput("SDEMarkTable_Batch"),
          downloadButton('SDEdownload', "Download table",
                         class = "btn-secondary"),
          htmlOutput("SDElegend_de_columns", container = h6)
        )
      )
    },
    
    
    ### Heatmaps of gene expression Panel ----
    tabPanel(
      "Heatmaps of gene expression",
      fluidPage(
        hr(),
        p(
          "Display a heatmap showing relative expression and proportion of cells expressing a gene or group of genes across all neurons. This function uses data from threshold 2. Color shows relative scaled expression for each gene across neuron types, and is not comparable between genes."
        ),
        warning_msg,
        textAreaInput(
          inputId = "HMgenelist",
          label = "Introduce a list of genes",
          value = "flp*\nWBGene00001447\nWBGene00001448,zig-4",
          width = "500px",
          height = "100px"
        ),
        fluidRow(
          column(width = 3,
                 selectInput(
                   inputId = "HMdataset",
                   label = "Choose dataset: Neurons (threshold 2), All cells (unfiltered)",
                   choices = c("Neurons only", "All cell types")
                 )
          ),
          column(width = 1, offset = -10,
                 checkboxInput(inputId = "HMreorder_rows",
                               label = "Reorder rows",
                               value = TRUE),
          )
        ),
        
        actionButton(
          "HMbutton_from_list",
          "Plot heatmap from list",
          class = "btn-info",
          icon = icon("hand-point-right")
          
          
        ),
        hr(),
        fluidRow(
          column(3,fileInput("HMfile_input", NULL,
                             accept = c(
                               "text/csv",
                               "text/comma-separated-values,text/plain",
                               "txt")
          )),
          #column(3,actionButton("resetFile", "Clear uploaded file")),
          column(3,   actionButton(
            "HMbutton_from_file",
            "Plot heatmap from file",
            class = "btn-secondary",
            icon = icon("hand-point-right")
          ))
        ),
        
        p(
          "You can identify circles by clicking on them."
        ),
        
        div(style="height:30px;width:800px;padding-left:10px;padding-right:10px;background-color:#ffffff;",fluidRow(verbatimTextOutput("vals", placeholder = TRUE))),
        #uiOutput("dynamic"),
        br(),
        br(),
        br(),
        plotOutput("heatmap", width = "100%", height = "auto", hover = "plot_hover"),
        
        hr(),
        downloadButton("downloadheatmap", "Download plot",
                     class = "btn-secondary"),
        downloadButton("downloadheatmap_svg", "Download plot as SVG",
                     class = "btn-secondary")
        
        
      )
    )
  )
)

