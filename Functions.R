
##### Modified Seurat functions
##### SCeNGEA APP April 2020
#options(repos = BiocManager::repositories())
#options("repos")

load("Dataset_6July_2021_noSeurat2.rda")

utr <- c("WBGene00023498","WBGene00023497","WBGene00004397","WBGene00006843",
         "WBGene00004010","WBGene00006789","WBGene00001135","WBGene00001079",
         "WBGene00001135","WBGene00006783","WBGene00000501","WBGene00006788",
         "WBGene00001555")




# Perform DE ----

# Dispatch
perform_de <- function(ident.1, ident.2, method, ...){

  cat("DE of ", ident.1," vs ",ident.2,"\n")
  # dispatch to proper test
  if(method == "Wilcoxon on single cells"){
    print("sc Wilcoxon")
    tableDEX <- perform_de_sc(ident.1 , ident.2, ...)
  } else if(method == "Pseudobulk: Wilcoxon"){
    print("sc Wilcoxon")
    tableDEX <- perform_de_pb_wilcoxon(ident.1 , ident.2, ...)
  } else if(method == "Pseudobulk: edgeR pairwise exact test"){
    print("sc Wilcoxon")
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
    dplyr::arrange(p_val_adj)
  
  tableDEX$p_val <-
    as.numeric(formatC(tableDEX$p_val, format = "e", digits = 3) %>% gsub(" ", "", .))
  tableDEX$p_val_adj <-
    as.numeric(formatC(
      tableDEX$p_val_adj,
      format = "e",
      digits = 3
    ) %>% gsub(" ", "", .))
  
  tableDEX
}

#~ single cell Wilcoxon ----
# note this is a rewrite of Seurat::FindMarkers that does only the minimum required here,
# and works directly with the content of the Seurat object( not a full Seurat object)
perform_de_sc <- function(ident.1 , ident.2, min.pct = 0.1, min.diff.pct = -Inf, logfc.threshold = 0.25){
  cells.1 <- allCells.metadata$Neuron.type %in% ident.1
  
  if(!is.null(ident.2)){
    cells.2 <- allCells.metadata$Neuron.type %in% ident.2
  } else{
    cells.2 <- !( allCells.metadata$Neuron.type %in% ident.1 )
  }
  
  
  expr.1 <- allCells.data[,cells.1]
  expr.2 <- allCells.data[,cells.2]
  
  features <- rownames(allCells.data)
  
  fc.results <- FoldChange(object = allCells.data,
                           slot = "data", 
                           cells.1 = which(cells.1), cells.2 = which(cells.2),
                           features = features,
                           mean.fxn = function(x) {
                             return(log(x = rowMeans(x = expm1(x = x)) + 1, 
                                        base = 2))
                           }, fc.name = "avg_log2FC",
                           pseudocount.use = 1, base = 2)
  
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
  
  
  data.use <- allCells.data[features, c(which(cells.1), which(cells.2)), drop = FALSE]
  j <- seq_len(sum(cells.1))
  
  
  p_val <- sapply(X = seq_along(features),
                  FUN = function(x) {
                    min(2 * min(limma::rankSumTestWithCorrelation(index = j, 
                                                                  statistics = data.use[x, ])), 1)
                  })
  
  de.results <- cbind(data.frame(gene = features,
                                 p_val = p_val),
                      fc.results[features, , drop = FALSE])
  de.results$p_val_adj = p.adjust(p = de.results$p_val, 
                                  method = "bonferroni", n = nrow(allCells.data))
  rownames(de.results) <- NULL
  
  de.results$avg_log2FC <-
    as.numeric(formatC(de.results$avg_log2FC, digits = 3) %>% gsub(" ", "", .))
  
  de.results
}



#~ pseudobulk Wilcoxon ----
perform_de_pb_wilcoxon <- function(ident.1, ident.2, ...){
  
  cols.group.1 <- which(startsWith(colnames(pseudobulk_matrix[,]), ident.1))
  cols.group.2 <- which(startsWith(colnames(pseudobulk_matrix[,]), ident.2))
  
  data.use <- pseudobulk_matrix[, c(cols.group.1,cols.group.2)]
  
  p_val <- sapply(X = seq_len(nrow(data.use)),
                  FUN = function(x) {
                    min(2 * min(limma::rankSumTestWithCorrelation(index = seq_along(cols.group.1), 
                                                                  statistics = data.use[x, ])), 1)
                  })
  
  p_val_adj <- p.adjust(p_val, method = "bonferroni")
  data.frame(gene = rownames(pseudobulk_matrix),
             p_val = p_val,
             p_val_adj = p_val_adj)
}



#~ pseudobulk edgeR ----
perform_de_pb_edger <- function(ident.1, ident.2, ...){
  
  et <- exactTest(edger_precomputed, pair = c(ident.1, ident.2))
  
  et$table |>
    tibble::rownames_to_column() |>
    dplyr::mutate(p_val_adj = p.adjust(PValue, method = "bonferroni")) |>
    dplyr::rename(gene = rowname,
                  p_val = PValue,
                  p_val_adj = p_val_adj)
}


# # Tests
# res3 <- perform_de("sc_wilcoxon","AVL", "AWC_OFF")
# res2 <- perform_de("pb_wilcoxon","AVL", "AWC_OFF")
# res1 <- perform_de("pb_edger","AVL", "AWC_OFF")
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


