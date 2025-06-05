#### CeNGEN-app 2018-2022
#### Report bugs at https://github.com/cengenproject/CengenApp/issues



server <- function(input, output, session) {
  
  
  
  observe({
    query_GET <- getQueryString()
    
    if(length(query_GET) > 0 &&
       identical(names(query_GET), "celltype") &&
         length(query_GET) == 1L &&
         query_GET %in% all_cell_types){
      
      message("Switching to Enriched tab")
      
      
      updateSelectInput(inputId = "Markers2",
                        selected = query_GET)
      
      updateTabsetPanel(inputId = "tabs",
                        selected = "Enriched Genes by cell type")
      
      
    }
    
  })
  
  
  
  
  
  ### Gene expression by cell type Panel ----
  
  #~ Tcell ----
  observeEvent(input$TCell, {
    
    cat("--> thresholded gene expression TCell\n")
    
    withProgress(message = "Obtaining information...", value = 0, {
      
      if( input$Tcell_cut == "All Cells Unfiltered" ) {
        
        load_as_needed("all.TPM.raw_th")
        th = all.TPM.raw_th
        
      } else {
        load_as_needed("ths")
        th = ths
      }
      
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
        
        t4 <- th |>
          dplyr::filter(threshold == input$Tcell_cut) |>
          dplyr::select(`Gene name` = gene_name,
                        `Gene ID` = id,
                        `Expression level` = input$Tcell_name) |>
          dplyr::arrange(desc( `Expression level` ),
                         `Gene name`) |>
          dplyr::mutate(`Expression level` = round(`Expression level`, digits = 3))
        
        t4_display <- t4 |>
          dplyr::select(-`Gene ID`)
        
        
        output$Tcell_name_table <-
          DT::renderDataTable({
            DT::datatable(
              t4_display,
              options = list(pageLength = 10, autoWidth = TRUE),
              rownames = FALSE,
              style = 'jQueryUI',
              class = 'cell-border stripe'
            ) %>% formatStyle(c(1:2), color = "black", background = 'white')
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
              write.csv(t4 , file, dec = ".", sep = "\t")
            }
          )
        
      } else {
        
        print("Tcell_name not in colnames(th)")
        
        output$Tcell_name_table  <- NULL
        
        output$Error1 <- isolate(renderText({
            shiny::validate(need(
              input$Tcell_name %in% colnames(th),
              message = paste0(
                "WARNING: If you want to query non-neuronal cells select All Cells Unfiltered"
              )
            ))
          }))
      }
      
    })
  })
  
  
  #~ Tgene ----
  observeEvent(input$TGene, {
    cat("--> thresholded gene expression TGene\n")
    
    
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
      
      load_as_needed("all.TPM.raw_th")
      th = all.TPM.raw_th
      
    } else {
      
      load_as_needed("ths")
      th = ths
      
    }
    
    
    withProgress(message = "Obtaining information...", value = 0, {
      
      if (var %in% dplyr::filter(th, threshold == input$Tgene_cut)$gene_name) {
        
        t3 <- th |>
          dplyr::filter(gene_name == var,
                        threshold == input$Tgene_cut) |>
          dplyr::select(-id, -threshold, -gene_name)
        
        t3 <- as.numeric(t3) |> setNames(colnames(t3))
        
        t3 <- sort(t3, decreasing = TRUE)
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
              is.null(unreliable_gene_names) ||
                (!dplyr::filter(gene_list, gene_name == var)$gene_name %in% unreliable_gene_names),
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
        print("No gene left")
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
          ) %>% formatStyle(c(1:2), color = "black", backgroundColor = 'white')
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
  
  #~ Tgene_batch ----
  observeEvent( list(input$Tgene_name_batch, input$Tgene_cut_batch) , {
    cat("--> thresholded gene expression Tgene_name_batch/Tgene_cut_batch\n")
    
    
    gns1 <- strsplit(as.character(input$Tgene_name_batch), "\n| |\\,|\t")
    gns1 <- as.data.frame(gns1)[,1]
    
    star <- grep("\\*", gns1)
    
    families1 <- c()
    for(i in star){ 
      families1 <- c(families1, grep( gsub("\\*", "", gns1[i]), gene_list$gene_name, value = TRUE) )
    }
    
    gns <- unique(c(families1, gns1, filter(gene_list, gene_id %in% gns1 | seqnames %in% gns1)$gene_name))
    
    
    if( input$Tgene_cut_batch == "All Cells Unfiltered" ) { 
      
      load_as_needed("all.TPM.raw_th")
      th = all.TPM.raw_th
      
    } else { 
      
      load_as_needed("ths")
      th = ths
      
    }
    
    if (length(which(gns %in% unique(th$gene_name))) > 0) {
      
      
      tb <- th |>
        dplyr::filter(gene_name %in% gns,
                      threshold == input$Tgene_cut_batch) |>
        dplyr::relocate(gene_name, .before = 1) |>
        dplyr::select(-id, -threshold)
      
      output$textb <- renderText({""})
      head(tb)
      req(input$Tgene_name_batch)
      
      
      output$TGeneBatch <-
        downloadHandler(
          filename = function() {
            paste0(
              "GenesExpressing-BATCH",
              "-thrs",
              input$Tgene_cut_batch,
              ".csv"
            )
          },
          content = function(file) {
            write.csv(tb, file)
          }
        )
      
    } else {
      print("No gene left")
      tb <- NULL
      output$textb <-
        renderText({
          "No gene found with this name in the dataset"
        })
    }
  })
  
  
  
  ### Find markers based on percentage of expression Panel ----
  
  observeEvent(input$Filter, {
    cat("--> find markers Filter\n")
    
    load_as_needed("pcttable")
    
    s1 <- unlist(strsplit(as.character(input$PCT_group1), split = ","))
    s1 <- gsub(" ", "", as.character(s1))
    s2 <- unlist(strsplit(as.character(input$PCT_group2), split = ","))
    s2 <- gsub(" ", "", as.character(s2))
    
    
    expressed_thr <- as.numeric( input$PCT_expressed_threshold )
    not_expressed_thr <- as.numeric( input$PCT_notExpressed_threshold )
    
    
    
    expr_in_grp1 <- pcttable |>
      dplyr::filter(id %in% s1,
                    pct.exp >= expressed_thr)
    
    genes_expr_in_grp1 <- names(which(table(expr_in_grp1$gene_name) == length(s1)))
    
    notexpr_in_grp2 <- pcttable |>
      dplyr::filter(id %in% s2,
                    pct.exp < not_expressed_thr)
    
    genes_notexpr_in_grp2 <- names(which(table(notexpr_in_grp2$gene_name) == length(s2)))
    
    
    genes_res <- intersect(genes_expr_in_grp1, genes_notexpr_in_grp2)
    
    
    res_table <- pcttable |>
      dplyr::filter(gene_name %in% genes_res) |>
      dplyr::select(`Gene Name` = gene_name,
                    `Gene ID` = gene_id) |>
      dplyr::distinct() |>
      as.data.frame()
    
    
    yes_display <- expr_in_grp1 |>
      dplyr::mutate(pct.exp = round(pct.exp, 2),
                    avg.exp = round(avg.exp, 2)) |>
      dplyr::arrange(desc(pct.exp))
    
    
    no_display <- notexpr_in_grp2 |>
      dplyr::mutate(pct.exp = round(pct.exp, 2),
                    avg.exp = round(avg.exp, 2)) |>
      dplyr::arrange(desc(pct.exp))
    
    yes_download <- dplyr::filter(expr_in_grp1, gene_name %in% genes_res)
    
    no_download <- dplyr::filter(notexpr_in_grp2, gene_name %in% genes_res)
    
    
    
    
    output$YesExpressed <-
      DT::renderDataTable({
        DT::datatable(
          yes_display,
          options = list(pageLength = 10, autoWidth = TRUE),
          rownames = FALSE,
          style = 'jQueryUI',
          class = 'cell-border stripe'
        ) %>% formatStyle(c(1:6), color = "black", backgroundColor = 'white')
      })
    output$NoExpressed <-
      DT::renderDataTable({
        DT::datatable(
          no_display,
          options = list(pageLength = 10, autoWidth = TRUE),
          rownames = FALSE,
          style = 'jQueryUI',
          class = 'cell-border stripe'
        ) %>% formatStyle(c(1:6), color = "black", backgroundColor = 'white')
      })
    output$Result <-
      DT::renderDataTable({
        DT::datatable(
          res_table,
          options = list(pageLength = 10, autoWidth = TRUE),
          rownames = FALSE,
          style = 'jQueryUI',
          class = 'cell-border stripe'
        ) %>% formatStyle(c(1:2), color = "black", backgroundColor = 'white')
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
          write.csv(rbind(yes_download, no_download), file, sep = "\t")
        }
      )
    
  })
  
  
  ### Enriched Genes by cell type Panel ----
  
  observeEvent(input$Markers, {
    cat("--> enriched genes Markers\n")
    
    load_as_needed("markers")
    
    print(input$top2)
    output$MarkTable <- DT::renderDataTable({
      DT::datatable(
        filter(markers, cluster == input$Markers, avg_log2FC > 0) %>%
          arrange(p_val_adj, desc(avg_log2FC)) %>%
          head(as.numeric(input$top)) ,
        style = 'jQueryUI',
        class = 'cell-border stripe',
        rownames = FALSE
      ) %>%
        formatStyle(c(1:8), color = "black", backgroundColor = 'white')
    })
    
    t1 <-
      filter(markers, cluster == input$Markers, avg_log2FC > 0) %>%
      arrange(p_val_adj, desc(avg_log2FC))
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
    cat("--> enriched genes Markers2\n")
    
    load_as_needed("markersAllcells")
    
    output$MarkTable2 <- DT::renderDataTable({
      DT::datatable(
        filter(markersAllcells, cluster == input$Markers2, avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% head(as.numeric(input$top2)),
        style = 'jQueryUI',
        class = 'cell-border stripe',
        rownames = FALSE
      ) %>% formatStyle(c(1:8), color = "black", backgroundColor = 'white')
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
  
  
  
  
  ### Find Differential Expression between Cell Types Panel ----
  observeEvent(input$DEbutton, {
    cat("--> testing DE of ", input$DEgroup1, " vs ", input$DEgroup2, "\n")
    
    b1 <- unlist(strsplit(input$DEgroup1, split = ","))
    b1 <- gsub(" ", "", as.character(b1))
    b2 <- unlist(strsplit(input$DEgroup2, split = ","))
    b2 <- gsub(" ", "", as.character(b2))
    
    b1_is_valid <- all(b1 %in% all_cell_types)
    b2_is_valid <- all(b2 %in% all_cell_types) ||
      "ALL" %in% b2 ||
      "NEURONS" %in% b2
    
    comparing_multiple_cell_types <- length(b1) > 1 ||
      length(b2) > 1||
      b2 == "ALL" ||
      b2 == "NEURONS"
    
    
    if (!b1_is_valid || !b2_is_valid) {
      
      output$MarkTable_Batch <- DT::renderDataTable({data.frame()})
      output$text_error_dex <- renderText({"One or more cell types introduced are not correct"})
      
    } else if(input$DEtest == "Pseudobulk: edgeR pairwise exact test" && comparing_multiple_cell_types){
      
      output$MarkTable_Batch <- DT::renderDataTable({data.frame()})
      output$text_error_dex <- renderText({"edgeR exact test can only compare pairs of cell types"})
      
    } else{
      
      
      if (any(b2 == "ALL")){
        print("testing vs ALL")
        
        b2 <- all_cell_types |> setdiff(b1)
        
      }
      
      if (any(b2 == "NEURONS")){
        print("Testing against all neurons.")
        
        load_as_needed("allCells.metadata")
        
        b2 <- all_neuron_types |> setdiff(b1) #does NOT contain "_stressed" neurons
      }
      
      tableDEX <-
        perform_de(
          ident.1 = b1,
          ident.2 = b2,
          method = input$DEtest
        )
      
      
      if(input$DEtest == "Pseudobulk: edgeR pairwise exact test" ||
         input$DEtest == "Pseudobulk: Wilcoxon"){
        
        load_as_needed("edger_precomputed")
        
        rows_b1 <- edger_precomputed$samples$group %in% b1
        rows_b2 <- edger_precomputed$samples$group %in% b2
        
        nb_sc_group_1 <- edger_precomputed$samples$nb_single_cells[rows_b1] |> sum()
        nb_sc_group_2 <- edger_precomputed$samples$nb_single_cells[rows_b2] |> sum()
        nb_rep_group_1 <- edger_precomputed$samples$nb_single_cells[rows_b1] |> length()
        nb_rep_group_2 <- edger_precomputed$samples$nb_single_cells[rows_b2] |> length()
        
        output$pseudobulk_metadata <- renderText({
          paste0("Comparing group 1 (",nb_sc_group_1," single cells in ",nb_rep_group_1,
                 " replicates) to group 2 (",nb_sc_group_2," single cells in ",
                 nb_rep_group_2," replicates)")
        }, sep = "<br>")
        
      } else{
        output$pseudobulk_metadata <- renderText({""}, sep = "<br>")
      }
      
      
      output$text_error_dex <- renderText({""})
      
      output$legend_de_columns <- renderText({
        
        if(input$DEtest == "Wilcoxon on single cells"){
          c("Testing differential expression between all single cells in the chosen clusters, using a Wilcoxon test. This test may display inflated power, as it considers each cell as an individual replicate.
          Before testing, the genes are filtered to only consider those displaying enough expression in one of the groups, 
          and a large enough fold change between groups.",
            "Columns:",
            "pct.1 and pct.2: percentage of single cells where the gene is detected in the first and second group.",
            "avg_logFC: expression change between group 1 and group 2 (as the log of the fold change of the means).",
            "p-val and FDR: nominal and adjusted (Benjamini-Hochberg) P-values of the test.")
          
        } else if(input$DEtest == "Pseudobulk: Wilcoxon"){
          c(
            "Testing differential expression between cell types across samples (pseudobulk), using a Wilcoxon test.
        Only genes which displayed high enough expression in enough cell types are considered.",
            "Columns:",
            "mean_1 and mean_2: mean expression across samples for group 1 and 2.",
            "log2FC: change in mean expression between group 1 and 2 (log fold change).",
            "p-val and FDR: nominal and adjusted (Benjamini-Hochberg) P-values of the test."
          )
        } else if(input$DEtest == "Pseudobulk: edgeR pairwise exact test"){
          c(
            "Testing differential expression between cell types across samples (pseudobulk), using edgeR's exact test.
        Only genes which displayed high enough expression in enough cell types are considered.
        Note that the exact test only allows pairwise comparisons, edgeR implements other tests (based on glm) which allow more complex comparisons but are impractical to run on a website.",
            "Columns:",
            "logFC: change in mean expression between group 1 and 2 (log fold change).",
            "logCPM: mean expression of the gene across all groups (log of Count Per Million reads)",
            "p-val and FDR: nominal and adjusted (Benjamini-Hochberg) P-values of the test."
          )
        }
      }, sep = "<br>")
      
      
      
      
      # Finalize output
      if (nrow(tableDEX) > 0) {
        
        output$MarkTable_Batch <- DT::renderDataTable({
          
          DT::datatable(
            tableDEX,
            style = 'jQueryUI',
            class = 'cell-border stripe',
            rownames = FALSE
          ) |>
            formatStyle(c(1:9), color = "black", backgroundColor = 'white') |>
            formatRound(columns = intersect(c("avg_logFC", "mean_1", "mean_2", "logCPM"),
                                            colnames(tableDEX)),
                                            digits = 1) |>
            formatSignif(columns =intersect(c("p_val", "FDR"),
                                            colnames(tableDEX)),
                                            digits = 2)
          
        })
        
        output$downloadDEX <-
          downloadHandler(
            filename = function() {
              paste(
                "DEgenes-",
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
        output$MarkTable_Batch <- DT::renderDataTable({NULL})
        output$text_error_dex <- renderText({"No feature passes the logfc threshold"})
      }
      
      
    }
  }, ignoreNULL = TRUE)
  
  
  
  
  
  
  ### Male vs Herm DE ----
  observeEvent(input$SDEbutton, {
    cat("--> sex-specific DE of ♂ ", input$SDEmale, " vs ⚥ ", input$SDEherm, "\n")
    
    b1 <- unlist(strsplit(input$SDEmale, split = ","))
    b1 <- gsub(" ", "", as.character(b1))
    b2 <- unlist(strsplit(input$SDEherm, split = ","))
    b2 <- gsub(" ", "", as.character(b2))
    
    b1_is_valid <- all(b1 %in% male_all_cell_types)
    b2_is_valid <- all(b2 %in% herm_all_cell_types)
    
    
    comparing_multiple_cell_types <- (length(b1) > 1) || (length(b2) > 1)
    
    
    if (!b1_is_valid || !b2_is_valid) {
      
      output$SDEMarkTable_Batch <- DT::renderDataTable({data.frame()})
      output$SDEtext_error_dex <- renderText({"One or more cell types introduced are not correct"})
      
    } else if(input$SDEtest == "Pseudobulk: edgeR pairwise exact test" && comparing_multiple_cell_types){
      
      output$SDEMarkTable_Batch <- DT::renderDataTable({data.frame()})
      output$SDEtext_error_dex <- renderText({"edgeR exact test can only compare pairs of cell types"})
      
    } else{
      
      
      tableDEX <-
        perform_de(
          ident.1 = b1,
          ident.2 = b2,
          method = input$SDEtest,
          subset_samples = c("male", "herm")
        )
      
      
      if(input$SDEtest == "Pseudobulk: edgeR pairwise exact test" ||
         input$SDEtest == "Pseudobulk: Wilcoxon"){
        
        load_as_needed("edger_precomputed")
        
        rows_b1 <- edger_precomputed$samples$group %in% b1
        rows_b2 <- edger_precomputed$samples$group %in% b2
        
        nb_sc_group_1 <- edger_precomputed$samples$nb_single_cells[rows_b1] |> sum()
        nb_sc_group_2 <- edger_precomputed$samples$nb_single_cells[rows_b2] |> sum()
        nb_rep_group_1 <- edger_precomputed$samples$nb_single_cells[rows_b1] |> length()
        nb_rep_group_2 <- edger_precomputed$samples$nb_single_cells[rows_b2] |> length()
        
        output$SDEpseudobulk_metadata <- renderText({
          paste0("Comparing group 1 (",nb_sc_group_1," single cells in ",nb_rep_group_1,
                 " replicates) to group 2 (",nb_sc_group_2," single cells in ",
                 nb_rep_group_2," replicates)")
        }, sep = "<br>")
        
      } else{
        output$SDEpseudobulk_metadata <- renderText({""}, sep = "<br>")
      }
      
      
      output$SDEtext_error_dex <- renderText({""})
      
      output$SDElegend_de_columns <- renderText({
        
        if(input$SDEtest == "Wilcoxon on single cells"){
          c("Testing differential expression between all single cells in the chosen clusters, using a Wilcoxon test.",
          "This test may display inflated power, as it considers each cell as an individual replicate.
          Before testing, the genes are filtered to only consider those displaying enough expression in one of the groups, 
          and a large enough fold change between groups.",
            "Columns:",
            "pct.1: percentage of single cells where the gene is detected in the male.",
          "pct.2: percentage of single cells where the gene is detected in the hermaphrodite",
            "avg_logFC: expression change between male and hermaphrodite (as the log of the fold change of the means). A positive number indicates higher expression in males.",
            "p-val and FDR: nominal and adjusted (Benjamini-Hochberg) P-values of the test.")
          
        } else if(input$SDEtest == "Pseudobulk: Wilcoxon"){
          c(
            "Testing differential expression between cell types across samples (pseudobulk), using a Wilcoxon test.
        Only genes which displayed high enough expression in enough cell types are considered.",
            "Columns:",
            "mean_1 and mean_2: mean expression across samples for males and hermaphrodites.",
            "log2FC: change in mean expression between males and hermaphrodites (log fold change). A positive number indicates higher expression in the male.",
            "p_val and FDR: nominal and adjusted (Benjamini-Hochberg) p-values of the test."
          )
        } else if(input$SDEtest == "Pseudobulk: edgeR pairwise exact test"){
          c(
            "Testing differential expression between cell types across samples (pseudobulk), using edgeR's exact test.
        Only genes which displayed high enough expression in enough cell types are considered.
        Note that the exact test only allows pairwise comparisons, edgeR implements other tests (based on glm) which allow more complex comparisons but are impractical to run on a website.",
            "Columns:",
            "logFC: change in mean expression between males and hermaphrodites (log fold change). A positive number indicates higher expression in males.",
            "logCPM: mean expression of the gene across all groups (log of Count Per Million reads)",
            "p_val and FDR: nominal and adjusted (Benjamini-Hochberg) p-values of the test."
          )
        }
      }, sep = "<br>")
      
      
      
      
      # Finalize output
      if (nrow(tableDEX) > 0) {
        
        output$SDEMarkTable_Batch <- DT::renderDataTable({
          
          DT::datatable(
            tableDEX,
            style = 'jQueryUI',
            class = 'cell-border stripe',
            rownames = FALSE
          ) |>
            formatStyle(c(1:9), color = "black", backgroundColor = 'white') |>
            formatRound(columns = intersect(c("avg_logFC", "mean_1", "mean_2", "logCPM"),
                                            colnames(tableDEX)),
                        digits = 1) |>
            formatSignif(columns =intersect(c("p_val", "FDR"),
                                            colnames(tableDEX)),
                         digits = 2)
          
        })
        
        output$SDEdownload <-
          downloadHandler(
            filename = function() {
              paste(
                "DEgenes-",
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
        
        output$SDEMarkTable_Batch <- DT::renderDataTable({NULL})
        
        if(! is.null(attr(tableDEX, "error_message")) ){
          
          output$SDEtext_error_dex <- renderText({attr(tableDEX, "error_message")})
        } else{
          
          output$SDEtext_error_dex <- renderText({"No feature passes the logfc threshold"})
        }
        
      }
      
      
    }
  }, ignoreNULL = TRUE)
  
  
  
  
  
  
  
  
  
  ### Heatmaps of gene expression Panel ----
  
  observeEvent(input$HMbutton_from_list, {
    
    cat("--> heatmap PlotHeatmapFromList\n")
    
    plot_heatmap(as.character(input$HMgenelist), input, output)
  })
  
  
  observeEvent(input$HMbutton_from_file, {
    
    cat("--> heatmap PlotHeatmapFromFile\n")
    
    inFile <- input$HMfile_input$datapath
    
    
    
    input_genelist <- read.table(inFile, header=FALSE)$V1
    
    plot_heatmap(input_genelist, input, output)
    
  })
  
}

