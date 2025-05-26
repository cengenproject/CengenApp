
##### Modified Seurat functions
##### SCeNGEA APP April 2020
#options(repos = BiocManager::repositories())
#options("repos")



# Loading data ----

# note the assignment in the global environment
load_as_needed <- function(dataset){
  if(!exists(dataset)){
    assign(dataset,
           qs::qread(paste0("data/2025-03-21/", dataset, ".qs")),
           envir = .GlobalEnv)
  }
}



# load global data
load_as_needed("gene_list")
load_as_needed("all_cell_types")
load_as_needed("all_neuron_types")







# Perform DE ----


#~ single cell Wilcoxon ----

# note these are rewrites of Seurat::FindMarkers and Seurat::FoldChange that do
# only the minimum required here, working directly with the content of
# the Seurat object (not a full Seurat object)

mean.fxn <- function(x) {
  return(log(x = rowMeans(x = expm1(x = x)) + 1, 
             base = 2))
}


# note this is a rewrite of Seurat::FoldChange that does only the minimum required here,
# and works directly with the content of the Seurat object( not a full Seurat object)
FoldChange <- function (cells.1, cells.2, features) {
  thresh.min <- 0
  pct.1 <- round(x = rowSums(x = allCells.data[features, cells.1, 
                                               drop = FALSE] > thresh.min)/length(x = cells.1), digits = 3)
  pct.2 <- round(x = rowSums(x = allCells.data[features, cells.2, 
                                               drop = FALSE] > thresh.min)/length(x = cells.2), digits = 3)
  data.1 <- mean.fxn(allCells.data[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(allCells.data[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c("avg_log2FC", "pct.1", "pct.2")
  return(fc.results)
}


perform_de_sc <- function(ident.1 , ident.2, min.pct = 0.1, min.diff.pct = -Inf, logfc.threshold = 0.25){
  
  
  load_as_needed("allCells.metadata")
  
  cells.1 <- allCells.metadata$Cell.type %in% ident.1
  
  cells.2 <- allCells.metadata$Cell.type %in% ident.2
  
  
  rm("allCells.metadata", envir = .GlobalEnv)
  
  load_as_needed("allCells.data")
  
  
  expr.1 <- allCells.data[,cells.1]
  expr.2 <- allCells.data[,cells.2]
  
  # get fold changes
  fc.results <- FoldChange(cells.1 = which(cells.1),
                           cells.2 = which(cells.2),
                           features = rownames(allCells.data))
  
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
perform_de_pb_wilcoxon <- function(ident.1, ident.2, ...){
  
  load_as_needed("pseudobulk_matrix")
  
  cols.group.1 <- which(startsWith(colnames(pseudobulk_matrix[,]), ident.1))
  cols.group.2 <- which(startsWith(colnames(pseudobulk_matrix[,]), ident.2))
  
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
             log2FC = log2FC,
             p_val = p_val,
             FDR = FDR) |>
    dplyr::arrange(FDR, p_val, desc(abs(log2FC))) |>
    mutate(p_val = signif(p_val, 2),
           FDR = signif(FDR, 2),
           avg_logFC = round(log2FC, 1),
           across(contains("mean"),
                  ~ round(.x, 1)))
}



#~ pseudobulk edgeR ----
perform_de_pb_edger <- function(ident.1, ident.2, ...){
  
  load_as_needed("edger_precomputed")
  
  et <- exactTest(edger_precomputed, pair = c(ident.2, ident.1))
  
  et$table |>
    tibble::rownames_to_column() |>
    dplyr::mutate(p_val_adj = p.adjust(PValue, method = "BH")) |>
    dplyr::rename(gene_id = rowname,
                  p_val = PValue,
                  FDR = p_val_adj) |>
    dplyr::arrange(FDR, p_val) |>
    dplyr::mutate(p_val = signif(p_val, 2),
                  FDR = signif(FDR, 2),
                  avg_logFC = round(logFC, 1))
}




#~ Dispatch ----
perform_de <- function(ident.1, ident.2, method, ...){
  
  cat("DE of ", ident.1," vs ",ident.2,"\n")
  # dispatch to proper test
  if(method == "Wilcoxon on single cells"){
    print("sc Wilcoxon")
    tableDEX <- perform_de_sc(ident.1 , ident.2, ...)
  } else if(method == "Pseudobulk: Wilcoxon"){
    print("pb Wilcoxon")
    tableDEX <- perform_de_pb_wilcoxon(ident.1 , ident.2, ...)
  } else if(method == "Pseudobulk: edgeR pairwise exact test"){
    print("pb edgeR")
    tableDEX <- perform_de_pb_edger(ident.1 , ident.2, ...)
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



