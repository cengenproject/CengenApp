
plot_heatmap <- function(input_gene_list, input, output){
  cat("--> heatmap PlotHeatmapFromList\n")
  
  
  
  #~~ Load dataset ----
  ds <- input$HMdataset
  
  if(ds == "Neurons only"){
    
    load_as_needed("med.scaled.long")
    
    heatmapdata = med.scaled.long |>
      select(gene_name, cell.type, scaled.expr, prop)
    
    
  } else {
    
    load_as_needed("TPM.raw.scaled.long")
    
    heatmapdata = TPM.raw.scaled.long |>
      select(gene_name, cell.type, scaled.expr, prop, tissue)
    
  }
  
  
  #~~ Process gene input list ----
  ss <- strsplit(input_gene_list, "\n| |\\,|\t")[[1]] |>
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
            paste(input_name_unknown, collapse = ","))
    message(length(input_known_but_no_expression)," genes known but no data: ",
            paste(input_known_but_no_expression, collapse = ","))
  }
  
  #~~ Determine gene order ----
  if ( length(ss_known) <= 1 || !input$HMreorder_rows ){
    ordered_ss_known <- ss_known
    
  } else {
    
    # get full expression matrix
    if(ds=="Neurons only"){
      
      load_as_needed("TPM.medium")
      expr_matrix <- TPM.medium
      
    } else {
      
      load_as_needed("all.TPM.raw")
      expr_matrix <- as(all.TPM.raw,"dgCMatrix")
      
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
    
    if(ds == "All cell types"){
      
      tissues_table <- heatmapdata[heatmapdata$gene_name == "nduo-6", c("cell.type", "tissue")]
      fake_heatmap_data <- merge(fake_heatmap_data, tissues_table, by = "cell.type")
    }
    
    heatmapdata <- rbind(heatmapdata, fake_heatmap_data)
  }
  
  
  #~~ Order cells and genes ----
  
  if(ds == "Neurons only"){
    
    heatmapdata$tissue <- "Neuron"
    
  }
  
  heatmapdata <- heatmapdata |>
    filter(gene_name %in% ordered_ss_known) |>
    mutate(gene_name = factor(gene_name,
                              levels = ordered_ss_known),
           special_order = case_when(
             cell.type %in% pharyngeal_neurons ~ 2,
             cell.type %in% male_neurons ~ 3,
             .default = 1
           ),
           tissue = factor(tissue,
                           levels = ordered_tissues)) |>
    arrange(tissue, special_order) |>
    mutate(cell.type = fct_inorder(cell.type)) |>
    select(- special_order) |>
    as.data.frame()
  
  
  
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
          strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          strip.clip = "off"
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
  
  
  output$vals <- renderPrint({
    hover <- input$plot_hover
    y <- nearPoints(heatmapdata, input$plot_hover)
    req(nrow(y) != 0)
    y
  })
}


