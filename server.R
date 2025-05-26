#### CeNGEN-app 2018-2022
#### Report bugs at https://github.com/cengenproject/CengenApp/issues



server <- function(input, output) {
  
  
  ### Gene expression by cell type Panel ----
  
  #~ Tcell ----
  observeEvent(input$TCell, {
    
    cat("--> thresholded gene expression TCell\n")
    
    withProgress(message = "Obtaining information...", value = 0, {
      
      if( input$Tcell_cut == "All Cells Unfiltered" ) {
        
        load_as_needed("L4.all.TPM.raw_th")
        th = L4.all.TPM.raw_th
        
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
      
      load_as_needed("L4.all.TPM.raw_th")
      th = L4.all.TPM.raw_th
      
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
              !dplyr::filter(gene_list, gene_name == var)$gene_id %in% unreliable_gene_ids,
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
      
      load_as_needed("L4.all.TPM.raw_th")
      th = L4.all.TPM.raw_th
      
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
  observeEvent(input$DEXButton, {
    cat("--> testing DE of ", input$batch1, " vs ", input$batch2, "\n")
    
    b1 <- unlist(strsplit(input$batch1, split = ","))
    b1 <- gsub(" ", "", as.character(b1))
    b2 <- unlist(strsplit(input$batch2, split = ","))
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
      
    } else if(input$test == "Pseudobulk: edgeR pairwise exact test" && comparing_multiple_cell_types){
      
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
          method = input$test
        )
      
      
      if(input$test == "Pseudobulk: edgeR pairwise exact test" ||
         input$test == "Pseudobulk: Wilcoxon"){
        
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
        
        if(input$test == "Wilcoxon on single cells"){
          c("Testing differential expression between all single cells in the chosen clusters, using a Wilcoxon test. This test may display inflated power, as it considers each cell as an individual replicate.
          Before testing, the genes are filtered to only consider those displaying enough expression in one of the groups, 
          and a large enough fold change between groups.",
            "Columns:",
            "pct.1 and pct.2: percentage of single cells where the gene is detected in the first and second group.",
            "avg_logFC: expression change between group 1 and group 2 (as the log of the fold change of the means).",
            "p-val and p_val_adj: nominal and adjusted (Bonferroni) P-values of the test.")
          
        } else if(input$test == "Pseudobulk: Wilcoxon"){
          c(
            "Testing differential expression between cell types across samples (pseudobulk), using a Wilcoxon test.
        Only genes which displayed high enough expression in enough cell types are considered.",
            "Columns:",
            "mean_1 and mean_2: mean expression across samples for group 1 and 2.",
            "log2FC: change in mean expression between group 1 and 2 (log fold change).",
            "p_val and p_val_adj: nominal and adjusted (Bonferroni) p-values of the test."
          )
        } else if(input$test == "Pseudobulk: edgeR pairwise exact test"){
          c(
            "Testing differential expression between cell types across samples (pseudobulk), using edgeR's exact test.
        Only genes which displayed high enough expression in enough cell types are considered.
        Note that the exact test only allows pairwise comparisons, edgeR implements other tests (based on glm) which allow more complex comparisons but are impractical to run on a website.",
            "Columns:",
            "logFC: change in mean expression between group 1 and 2 (log fold change).",
            "logCPM: mean expression of the gene across all groups (log of Count Per Million reads)",
            "p_val and p_val_adj: nominal and adjusted (Bonferroni) p-values of the test."
          )
        }
      }, sep = "<br>")
      
      
      
      
      # Finalize output
      if (nrow(tableDEX) > 0) {
        
        output$MarkTable_Batch <- DT::renderDataTable({
          DT::datatable(
            tableDEX |> head( as.numeric(input$topM2) ),
            options = list( pageLength = as.numeric(input$topM2) ),
            style = 'jQueryUI',
            class = 'cell-border stripe',
            rownames = FALSE
          ) |> formatStyle(c(1:8), color = "black", backgroundColor = 'white') |>
            formatRound(columns = c('avg_logFC'), digits = 1) |>
            formatSignif(columns = c('p_val', 'p_val_adj'), digits = 2)
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
        output$MarkTable_Batch <- DT::renderDataTable({NULL})
        output$text_error_dex <- renderText({"No feature passes the logfc threshold"})
      }
      
      
    }
  }, ignoreNULL = TRUE)
  
  
  
  
  
  
  
  
  
  ### Heatmaps of gene expression Panel ----
  
  #~ From list ----
  observeEvent(input$PlotHeatmapFromList, {
    
    cat("--> heatmap PlotHeatmapFromList\n")
    
    
    
    #~~ Load dataset ----
    ds <- input$dataset_heatmap
    
    if(ds == "Neurons only"){
      
      load_as_needed("med.scaled.long")
      
      heatmapdata=med.scaled.long
      
      
    } else {
      
      load_as_needed("L4.TPM.raw.scaled.long")
      
      heatmapdata=L4.TPM.raw.scaled.long
      
    }
    
    
    #~~ Process gene input list ----
    ss <- strsplit(as.character(input$genelist), "\n| |\\,|\t")[[1]] |>
      setdiff("") |>
      unique()
    
    
    # expand stars
    input_has_star <- grepl("\\*", ss)
    
    expanded_stars <- lapply(which(input_has_star),
                             \(.i){
                               pattern <- paste0("^", gsub("\\*", ".+", ss[[.i]]), "$")
                               grep(pattern, gene_list$gene_name, value = TRUE)
                             }) |>
      unlist()
    
    ss <- ss[!input_has_star]
    ss <- union(ss, expanded_stars)
    
    
    # convert to gene_names
    is_a_gene_id <- !(ss %in% gene_list$gene_name) &
      ss %in% gene_list$gene_id
    
    is_a_seq_name <- !(ss %in% gene_list$gene_name) &
      ss %in% gene_list$seqnames
    
    ss[is_a_gene_id] <- vlookup(ss[is_a_gene_id], gene_list, result_column = 2, lookup_column = 1)
    ss[is_a_seq_name] <- vlookup(ss[is_a_seq_name], gene_list, result_column = 2, lookup_column = 3)
    
    
    
    input_name_unknown <- ss[ !(ss %in% gene_list$gene_name) ]
    input_has_no_expression_data <- ss[!ss %in% heatmapdata$gene_name]
    
    input_known_but_no_expression <- setdiff(input_has_no_expression_data, input_name_unknown)
    
    ss_known <- ss |> setdiff(input_has_no_expression_data)
    
    message("Heatmap from list, ", length(ss_known)," genes: ", paste(ss_known, collapse = ","))
    if(length(input_has_no_expression_data) > 0){
      message(length(input_name_unknown)," genes unknown: ",
              input_name_unknown)
      message(length(input_known_but_no_expression)," genes known but no data: ",
              input_known_but_no_expression)
    }
    
    #~~ Reorder ----
    if ( length(ss_known) <= 1 || !input$reorder_rows ){
      ordered_ss_known <- ss_known
      
    } else {
      
      # get full expression matrix
      if(ds=="Neurons only"){
        
        load_as_needed("L4.TPM.medium")
        expr_matrix <- L4.TPM.medium
        
      } else {
        
        load_as_needed("L4.all.TPM.raw")
        expr_matrix <- as(L4.all.TPM.raw,"dgCMatrix")
        
      }
      
      ss_known_gene_id <- vlookup(ss_known, gene_list, result_column = 1, lookup_column = 2)
      
      order_for_ss_ids <- expr_matrix[ss_known_gene_id, , drop = FALSE] |>
        t() |>
        scale() |>
        t() |>
        dist() |>
        hclust() |>
        purrr::pluck("order")
      
      ordered_ss_known <- ss_known[order_for_ss_ids]
      
    }
    
    # final plotting order
    ordered_ss_known <- c(
      rev(ordered_ss_known),
      input_known_but_no_expression
    )
    
    
    
    #~~ Prepare dataset for plotting ----
    
    
    if(length(input_known_but_no_expression) > 0){
      
      cell_types <- heatmapdata$cell.type[heatmapdata$gene_name == "nduo-6"]
      
      fake_heatmap_data <- expand.grid(
        gene_name=input_known_but_no_expression,
        cell.type= cell_types,
        scaled.expr=0,
        prop=0
      )
      
      if(ds == "Neurons only"){
        
        fake_heatmap_data$Modality <- "NA"
        
      } else{
        
        tissues_table <- heatmapdata[heatmapdata$gene_name == "nduo-6", c("cell.type", "tissue")]
        fake_heatmap_data <- merge(fake_heatmap_data, tissues_table, by = "cell.type")
      }
      
      heatmapdata <- rbind(heatmapdata, fake_heatmap_data)
    }
    
    heatmapdata <- heatmapdata[heatmapdata$gene_name %in% ordered_ss_known, ]
    
    
    heatmapdata$gene_name <- factor(heatmapdata$gene_name,
                                    levels = ordered_ss_known)
    
    
    
    
    #~~ Plot ----
    
    g_base <- ggplot(heatmapdata,
                     aes(y = gene_name, x = cell.type)) +
      geom_point(aes(color = scaled.expr, size = prop)) +
      scale_color_gradientn("Scaled TPM", colors = c("orange", "maroon", "navy")) +
      scale_size_continuous(name = "Proportion", limit = c(0.5, 100), range = c(0,5)) +
      theme(
        panel.background = element_blank()
      )
    
    
    
    if (ds != "Neurons only"){
      
      g <- g_base +
        labs(y = "Gene", x = "Tissue") +
        theme(
          axis.text.y.left = element_text(size = 12, color = "black"),
          legend.key = element_blank(),
          legend.text = element_text(size = 14, color = "black"),
          legend.title = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          # grid
          panel.grid = element_line(linewidth = 0.2, color = "grey85")
        )
      
      
      if(length(ordered_ss_known) > 1){
        g <- g +
          facet_grid(~tissue, scales = "free_x", space = "free_x", switch = "x") +
          theme(
            strip.placement = "outside", 
            strip.background.x = element_blank(),
            axis.text.x.bottom = element_blank(),
            strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)
          ) 
      }
      
      pg <- ggplotGrob(g)
      
      
      for(i in which(grepl("strip-b", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
      }
      
      fnh <-
        function() {
          withProgress(message = "Generating heatmap Plot...", value = 0, {
            pg      })
        }
      
      
    } else {
      
      g <- g_base +
        labs(y = "Gene", x= "Neuron") +
        theme(
          axis.text.x.bottom = element_text(angle = 90, vjust= 0.5,
                                            size = 7, color = "black",
                                            hjust = 1),
          axis.text.y.left = element_text(size = 7, color = "black"),
          panel.grid = element_line(linewidth = 0.5, color = "grey85")
        )
      
      
      fnh <-
        function() {
          withProgress(message = "Generating heatmap Plot...", value = 0, { g })
        }
    }
    
    l <- length(unique(g$data$gene_name))
    
    output$downloadheatmap <-
      downloadHandler(
        filename = function() {
          filename = paste("Heatmap.png", sep = "")
        },
        content = function(file) {
          ggsave(
            fnh(),
            file = file,
            height = 1500 + 20*l,
            width = 5000,
            units = "px" ,
            limitsize = FALSE,
            device = "png"
          )
        }
      )
    
    output$downloadheatmap_svg <-
      downloadHandler(
        filename = function() {
          filename = paste("Heatmap.svg", sep = "")
        },
        content = function(file) {
          ggsave(
            fnh(),
            file = file,
            height = 1500 + 20*l,
            width = 5000,
            units = "px" ,
            limitsize = FALSE,
            device = "svg"
          )
        }
      )
    
    output$heatmap <- renderPlot(g, height = 250 + 10*l)
    
    output$dynamic <- renderUI({
      #req(input$plot_hover)
      verbatimTextOutput("vals", placeholder = TRUE)
    })
    
    output$vals <- renderPrint({
      hover <- input$plot_hover 
      #print(input$plot_hover) # list
      y <- nearPoints(heatmapdata, input$plot_hover)
      req(nrow(y) != 0)
      y
    })
  })
  
  #~ From file ----
  
  observeEvent(input$PlotHeatmapFromFile, {
    cat("--> heatmap PlotHeatmapFromFile\n")
    
    load_as_needed("med.scaled.long")
    load_as_needed("L4.TPM.raw.scaled.long")
    
    ds <- input$dataset_heatmap
    inFile <- input$file1
    ss<-read.table(inFile$datapath, header=FALSE)$V1
    
    # ss <- strsplit(as.character(input$genelist), "\n| |\\,|\t")
    # ss <- as.data.frame(ss)[,1]
    star <- grep("\\*", ss)
    
    families <- c()
    for(i in star){ 
      families <- c(families, grep( gsub("\\*", "", ss[i]), gene_list$gene_name, value = TRUE) )
    }
    
    ss <- unique(c(families, ss, filter(gene_list, gene_id %in% ss | seqnames %in% ss)$gene_name))
    
    mis <- ss[ss %in% c(gene_list$gene_id, gene_list$gene_name, gene_list$seqnames) & !ss %in% med.scaled.long$gene_name & ss %in% gene_list$gene_name]
    mis_all <- ss[ss %in% c(gene_list$gene_id, gene_list$gene_name, gene_list$seqnames) & !ss %in% L4.TPM.raw.scaled.long$gene_name & ss %in% gene_list$gene_name]
    
    if(ds=="Neurons only"){
      
      load_as_needed("L4.TPM.medium")
      load_as_needed("ths")
      load_as_needed("L4.all.TPM.raw")
      
      L4.TPM=L4.TPM.medium
      heatmapdata=med.scaled.long
      cc = colnames(ths)[-c(1,130,131)]
      missing = mis
    } else {
      
      load_as_needed("L4.all.TPM.raw")
      
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
    
    # order
    if ( nrow(flp.expr) > 1 && input$reorder_rows ) {
      flp.neuron.order <- pheatmap(flp.expr, scale = "row")
      flp.neuron.order <- flp.neuron.order[["tree_row"]]$order
      flp.neuron.order <- flp.ids[flp.neuron.order]
      
    } else if(nrow(flp.expr) > 1 && ! input$reorder_rows){
      flp.neuron.order <- as.character(vlookup(unique(ss), gene_list, result_column = 1, lookup_column = 2))
      
    } else{
      flp.neuron.order <- flp.ids[1]
    }
    
    flp.neuron.order <- as.character(vlookup(flp.neuron.order, gene_list))
    flp.neuron.scaled$gene_name <- factor(flp.neuron.scaled$gene_name, levels = c(rev(flp.neuron.order), missing))
    #flp.neuron.scaled$gene_name <- fct_rev(flp.neuron.scaled$gene_name)
    
    for( i in missing ){
      dff <- data.frame(gene_name=i, cell.type= cc, scaled.expr=0, prop=0, Modality="NA")
      if(ds!="Neurons only"){
        colnames(dff)[5]<-"tissue"
        dff$tissue <- filter(L4.TPM.raw.scaled.long, gene_name=="nduo-6")$tissue
      }
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
            pg      })
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
    
    l <- length(unique(g$data$gene_name))
    
    output$downloadheatmap <-
      downloadHandler(
        filename = function() {
          filename = paste("Heatmap.png", sep = "")
        },
        content = function(file) {
          ggsave(
            fnh(),
            file = file,
            height = 1500 + 20*l,
            width = 5000,
            units = "px",
            limitsize = FALSE,
            device = "png"
          )
        }
      )
    
    output$downloadheatmap_svg <-
      downloadHandler(
        filename = function() {
          filename = paste("Heatmap.svg", sep = "")
        },
        content = function(file) {
          ggsave(
            fnh(),
            file = file,
            height = 1500 + 20*l,
            width = 5000,
            units = "px" ,
            limitsize = FALSE,
            device = "svg"
          )
        }
      )
    
    output$heatmap <- renderPlot(g, height = 250 + 10*l)
    
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