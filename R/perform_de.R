

#~ single cell Wilcoxon ----

# note these are rewrites of Seurat::FindMarkers and Seurat::FoldChange that do
# only the minimum required here, working directly with the content of
# the Seurat object (not a full Seurat object)

mean.fxn <- function(x) {
  log(rowMeans(expm1(x)) + 1,
      base = 2)
}


# note this is a rewrite of Seurat::FoldChange that does only the minimum required here,
# and works directly with the content of the Seurat object( not a full Seurat object)
FoldChange <- function (cells.1, cells.2, cnt_mat) {
  
  thresh.min <- 0
  
  thres_mat_1 <- cnt_mat[, cells.1, drop = FALSE] > thresh.min
  thres_mat_2 <- cnt_mat[, cells.2, drop = FALSE] > thresh.min
  
  pct.1 <- round(x = rowSums(thres_mat_1)/length(cells.1), digits = 3)
  pct.2 <- round(x = rowSums(thres_mat_2)/length(cells.2), digits = 3)
  
  data.1 <- mean.fxn(thres_mat_1)
  data.2 <- mean.fxn(thres_mat_2)
  
  fc <- (data.1 - data.2)
  
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c("avg_log2FC", "pct.1", "pct.2")
  
  return(fc.results)
}


perform_de_sc <- function(ident.1 , ident.2, min.pct = 0.1, min.diff.pct = -Inf,
                          logfc.threshold = 0.25, subset_samples = NULL){
  
  
  if(is.null(subset_samples)){
    
    load_as_needed("allCells.metadata")
    
    cells.1 <- allCells.metadata$Cell.type %in% ident.1
    cells.2 <- allCells.metadata$Cell.type %in% ident.2
    
    rm("allCells.metadata", envir = .GlobalEnv)
    
    # do not rely on `load_as_needed()` since we don't know what the previous state was
    allCells.data <- qs::qread(file.path(data_dir, "allCells.data.qs"))
    
  } else{
    
    stopifnot(length(subset_samples) == 2)
    
    load_as_needed("comb_allCells.metadata")
    
    cells.1 <- ( comb_allCells.metadata$Cell.type %in% ident.1 ) &
      ( comb_allCells.metadata$sample_set == subset_samples[[1]] )
    
    cells.2 <- ( comb_allCells.metadata$Cell.type %in% ident.2 ) &
      ( comb_allCells.metadata$sample_set == subset_samples[[2]] )
    
    rm("comb_allCells.metadata", envir = .GlobalEnv)
    
    
    # do not rely on `load_as_needed()` since we don't know what the previous state was
    allCells.data <- qs::qread(file.path(data_dir, "comb_allCells.data.qs"))
    
  }
  
  
  expr.1 <- allCells.data[,cells.1]
  expr.2 <- allCells.data[,cells.2]
  
  # get fold changes
  fc.results <- FoldChange(cells.1 = which(cells.1),
                           cells.2 = which(cells.2),
                           cnt_mat = allCells.data)
  
  
  
  # filter features
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- rownames(fc.results)[alpha.min >= min.pct &
                                     alpha.diff >= min.diff.pct]
  
  if (length(features) == 0) {
    warning("No features pass min threshold; returning empty data.frame")
  }
  
  features.diff <- rownames(fc.results)[abs(fc.results[["avg_log2FC"]]) >= logfc.threshold]
  
  features <- intersect(features, features.diff)
  if (length(features) == 0) {
    warning("No features pass logfc.threshold threshold; returning empty data.frame")
  }
  
  
  # Wilcoxon test
  data.use <- allCells.data[features, c(which(cells.1), which(cells.2)), drop = FALSE]
  j <- seq_len(sum(cells.1))
  
  p_val <- sapply(X = seq_along(features),
                  FUN = function(x) {
                    min(2 * min(limma::rankSumTestWithCorrelation(index = j, 
                                                                  statistics = data.use[x, ])), 1)
                  })
  
  
  p_val_adj = p.adjust(p_val, method = "BH", n = nrow(allCells.data))
  
  data.frame(gene_id = features,
             pct.1 = fc.results[features,]$pct.1,
             pct.2 = fc.results[features,]$pct.2,
             avg_logFC = fc.results[features,]$avg_log2FC,
             p_val = p_val,
             FDR = p_val_adj) |>
    dplyr::arrange(p_val_adj, p_val, desc(abs(avg_logFC)))
}



#~ pseudobulk Wilcoxon ----
perform_de_pb_wilcoxon <- function(ident.1, ident.2, subset_samples){
  
  if(is.null(subset_samples)){
    
    pseudobulk_matrix <- qs::qread(file.path(data_dir, "pseudobulk_matrix.qs"))
    
    cols.group.1 <- which(startsWith(colnames(pseudobulk_matrix), ident.1))
    cols.group.2 <- which(startsWith(colnames(pseudobulk_matrix), ident.2))
    
    
  } else{
    
    pseudobulk_matrix <- qs::qread(file.path(data_dir, "comb_pseudobulk_matrix.qs"))
    
    sample_set <- attr(pseudobulk_matrix, "sample_set")
    
    cols.group.1 <- which(startsWith(colnames(pseudobulk_matrix), ident.1) &
                            sample_set == subset_samples[[1]] )
    cols.group.2 <- which(startsWith(colnames(pseudobulk_matrix), ident.2) &
                            sample_set == subset_samples[[2]] )
    
  }
  
  
  
  idents <- c(ident.1, ident.2)
  missing_idents <- c(length(cols.group.1) == 0, length(cols.group.2) == 0)
  
  if(any(missing_idents)) {
    
    error_msg <- paste0("Error: cell type(s) not available: ", 
                        paste(idents[missing_idents], collapse = " and "))
    
    res <- data.frame(gene_id = character(0))
    attr(res, "error_message") <- error_msg
    
    return(res)
  }
  
  
  data.use <- pseudobulk_matrix[, c(cols.group.1,cols.group.2)]
  
  # Get fold change
  mean_1 <- rowMeans(pseudobulk_matrix[, cols.group.1])
  mean_2 <- rowMeans(pseudobulk_matrix[, cols.group.2])
  log2FC <- log2(mean_1 + 1) - log2(mean_2 + 1)
  
  
  p_val <- sapply(X = seq_len(nrow(data.use)),
                  FUN = function(x) {
                    min(2 * min(limma::rankSumTestWithCorrelation(index = seq_along(cols.group.1), 
                                                                  statistics = data.use[x, ])), 1)
                  })
  
  FDR <- p.adjust(p_val, method = "BH")
  
  data.frame(gene_id = rownames(pseudobulk_matrix),
             mean_1 = mean_1,
             mean_2 = mean_2,
             avg_logFC = log2FC,
             p_val = p_val,
             FDR = FDR) |>
    dplyr::arrange(FDR, p_val, desc(abs(avg_logFC)))
}



#~ pseudobulk edgeR ----
perform_de_pb_edger <- function(ident.1, ident.2, subset_samples = subset_samples){
  
  
  if(is.null(subset_samples)){
    
    edger_precomputed <- qs::qread(file.path(data_dir, "edger_precomputed.qs"))
    
    
  } else{
    
    edger_precomputed <- qs::qread(file.path(data_dir, "comb_edger_precomputed.qs"))
    
    ident.1 <- paste0(ident.1, "_", str_to_sentence(subset_samples[[1]]))
    
    ident.2 <- paste0(ident.2, "_", str_to_sentence(subset_samples[[2]]))
    
  }
  
  
  idents <- c(ident.1, ident.2)
  missing_idents <- !(idents %in% edger_precomputed$samples$group)
  
  if(any(missing_idents)) {
    
    error_msg <- paste0("Error: cell type(s) not available: ", 
                        paste(idents[missing_idents], collapse = " and "))
    
    res <- data.frame(gene_id = character(0))
    attr(res, "error_message") <- error_msg
    
    return(res)
  }
  
  
  et <- exactTest(edger_precomputed, pair = c(ident.2, ident.1))
  
  et$table |>
    tibble::rownames_to_column() |>
    dplyr::mutate(p_val_adj = p.adjust(PValue, method = "BH")) |>
    dplyr::rename(gene_id = rowname,
                  p_val = PValue,
                  FDR = p_val_adj,
                  avg_logFC = logFC) |>
    dplyr::arrange(FDR, p_val)
}




#~ Dispatch ----
perform_de <- function(ident.1, ident.2, method, subset_samples = NULL){
  
  cat("DE of ", ident.1," vs ",ident.2,"\n")
  # dispatch to proper test
  if(method == "Wilcoxon on single cells"){
    
    print("sc Wilcoxon")
    tableDEX <- perform_de_sc(ident.1 , ident.2, subset_samples = subset_samples)
    
  } else if(method == "Pseudobulk: Wilcoxon"){
    
    print("pb Wilcoxon")
    tableDEX <- perform_de_pb_wilcoxon(ident.1 , ident.2, subset_samples = subset_samples)
    
  } else if(method == "Pseudobulk: edgeR pairwise exact test"){
    
    print("pb edgeR")
    tableDEX <- perform_de_pb_edger(ident.1 , ident.2, subset_samples = subset_samples)
    
  } else{
    
    print("Test not recognized: ", method)
    stop("Test not recognized: ", method)
  }
  
  # finish
  
  left_join(tableDEX, gene_list, by = c("gene_id")) |>
    dplyr::relocate(gene_name, gene_id, -seqnames)
}


# # Tests
# res3 <- perform_de("AVL", "AWC_OFF", "Wilcoxon on single cells")
# res2 <- perform_de("AVL", "AWC_OFF", "Pseudobulk: Wilcoxon")
# res1 <- perform_de("AVL", "AWC_OFF", "Pseudobulk: edgeR pairwise exact test")
# 
# head(res1)
# head(res2)
# head(res3)
# 
# list(wilcox = res2$gene[res2$FDR < .05],
#      edgeR = res1$gene[res1$FDR < .05],
#      sc = res3$gene[res3$p_val_adj < .05]) |>
#   eulerr::euler() |>
#   plot()
# 
# left_join(res3 |>
#                   select(gene, sc_wilcox = avg_logFC),
#                 res1 |>
#                   select(gene, pb_edger = logFC),
#                 by = "gene") |>
#   ggplot() +
#   geom_point(aes(x = sc_wilcox, y = pb_edger))
# 
# left_join(res3 |>
#             select(gene, sc_wilcox = avg_logFC),
#           res2 |>
#             select(gene, pb_wilcox = log2FC),
#           by = "gene") |>
#   ggplot() +
#   geom_point(aes(x = sc_wilcox, y = pb_wilcox))
# left_join(res1 |>
#             select(gene, pb_edger = logFC),
#           res2 |>
#             select(gene, pb_wilcox = log2FC),
#           by = "gene") |>
#   ggplot() +
#   geom_point(aes(x = pb_edger, y = pb_wilcox))


