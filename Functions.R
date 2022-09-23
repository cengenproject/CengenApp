
##### Modified Seurat functions
##### SCeNGEA APP April 2020
#options(repos = BiocManager::repositories())
#options("repos")



# Loading data ----

# note the assignment in the global environment
load_as_needed <- function(dataset){
  if(!exists(dataset)){
    assign(dataset,
           qs::qread(paste0("data/", dataset, ".qs")),
           envir = .GlobalEnv)
  }
}

# load("Dataset_6July_2021_noSeurat2.rda")

# load global data
load_as_needed("gene_list")
load_as_needed("all_cell_types")
load_as_needed("all_neuron_types")


utr <- c("WBGene00023498","WBGene00023497","WBGene00004397","WBGene00006843",
         "WBGene00004010","WBGene00006789","WBGene00001135","WBGene00001079",
         "WBGene00001135","WBGene00006783","WBGene00000501","WBGene00006788",
         "WBGene00001555")





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
  
  load_as_needed("allCells.data")
  load_as_needed("allCells.metadata")
  
  cells.1 <- allCells.metadata$Neuron.type %in% ident.1
  
  cells.2 <- allCells.metadata$Neuron.type %in% ident.2
  
  
  
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
  
  
  p_val_adj = p.adjust(p_val, method = "bonferroni", n = nrow(allCells.data))
 
  data.frame(gene = features,
             pct.1 = fc.results$pct.1[features],
             pct.2 = fc.results$pct.2[features],
             avg_logFC = fc.results$avg_log2FC[features],
             p_val = p_val,
             p_val_adj = p_val_adj)
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
  
  p_val_adj <- p.adjust(p_val, method = "bonferroni")
  data.frame(gene = rownames(pseudobulk_matrix),
             mean_1 = mean_1,
             mean_2 = mean_2,
             log2FC = log2FC,
             p_val = p_val,
             p_val_adj = p_val_adj)
}



#~ pseudobulk edgeR ----
perform_de_pb_edger <- function(ident.1, ident.2, ...){
  
  load_as_needed("edger_precomputed")
  
  et <- exactTest(edger_precomputed, pair = c(ident.1, ident.2))
  
  et$table |>
    tibble::rownames_to_column() |>
    dplyr::mutate(p_val_adj = p.adjust(PValue, method = "bonferroni")) |>
    dplyr::rename(gene = rowname,
                  p_val = PValue,
                  p_val_adj = p_val_adj)
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
  
  tableDEX <-
    merge(
      tableDEX,
      gene_list,
      by.x = "gene",
      by.y = "gene_id",
      all.x = TRUE
    ) |>
    dplyr::arrange(p_val_adj) |>
    mutate(across(starts_with("p_val"),
                  ~ signif(.x, 2)),
           across(contains("log"),
                  ~ round(.x, 1)),
           across(contains("mean"),
                  ~ round(.x, 1)))
  
  
  tableDEX
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
# list(wilcox = res2$gene[res2$p_val_adj < .05],
#      edgeR = res1$gene[res1$p_val_adj < .05],
#      sc = res3$gene[res3$p_val_adj < .05]) |>
#   eulerr::euler() |>
#   plot()


