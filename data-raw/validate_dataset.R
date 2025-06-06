# if any of this fails, there will be harder to diagnose problems in the app.

## Inits ----

# ordered_tissues <- local({
#   xx <- parse("global.R")
#   eval(xx[sapply(xx, \(expr) expr[[2]] == "ordered_tissues")])
#   ordered_tissues
# })
# 
# validate_dataset("data/2025-03-21_L4/")
# validate_dataset("data/2025-03-17_L1/")
# validate_dataset("data/052225_herm/")
# validate_dataset("data/052225_male/male/")




validate_dataset <- function(data_dir){
  
  # File list ----
  expected_files <- c(
    "all_cell_types.qs", "all_neuron_types.qs", "gene_list.qs",
    "all.TPM.raw_th.qs", "ths.qs",
    "all.TPM.raw.qs", "TPM.medium.qs",
    "TPM.raw.scaled.long.qs", "med.scaled.long.qs",
    "pcttable.qs", "markers.qs", "markersAllcells.qs",
    "allCells.data.qs", "allCells.metadata.qs", "edger_precomputed.qs", "pseudobulk_matrix.qs"
  )
  
  stopifnot(all(
    expected_files %in% list.files(data_dir)
  ))
  
  
  
  all_cell_types <- qs::qread(file.path(data_dir, "all_cell_types.qs"))
  all_neuron_types <- qs::qread(file.path(data_dir, "all_neuron_types.qs"))
  gene_list <- qs::qread(file.path(data_dir, "gene_list.qs"))
  
  
  
  
  # thresholded ----
  
  #~ all.TPM.raw_th ----
  all.TPM.raw_th <- qs::qread(file.path(data_dir, "all.TPM.raw_th.qs"))
  
  
  stopifnot(identical(
    all_cell_types |> c("id", "gene_name", "threshold"),
    colnames(all.TPM.raw_th)
  ))
  
  
  stopifnot(all(
    all.TPM.raw_th$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    all.TPM.raw_th$id %in% gene_list$gene_id
  ))
  
  
  
  #~ ths ----
  ths <- qs::qread(file.path(data_dir, "ths.qs"))
  
  stopifnot(identical(
    all_neuron_types |> c("id", "gene_name", "threshold") |> sort(),
    colnames(ths) |> sort()
  ))
  
  
  stopifnot(all(
    ths$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    ths$id %in% gene_list$gene_id
  ))
  
  
  # Heatmap ----
  
  #~ med.scaled.long ----
  
  med.scaled.long <- qs::qread(file.path(data_dir, "med.scaled.long.qs"))
  
  stopifnot(all(
    c("gene_name", "cell.type", "scaled.expr", "prop") %in% colnames(med.scaled.long)
  ))
  
  stopifnot(all(
    med.scaled.long$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    med.scaled.long$cell.type %in% all_neuron_types
  ))
  
  
  #~ TPM.raw.scaled.long ----
  
  TPM.raw.scaled.long <- qs::qread(file.path(data_dir, "TPM.raw.scaled.long.qs"))
  
  stopifnot(all(
    c("gene_name", "cell.type", "scaled.expr", "prop", "tissue") %in% colnames(TPM.raw.scaled.long)
  ))
  
  stopifnot(all(
    TPM.raw.scaled.long$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    TPM.raw.scaled.long$cell.type %in% all_cell_types
  ))
  
  stopifnot(all(
    unique(TPM.raw.scaled.long$tissue) %in% ordered_tissues
  ))
  
  
  # Markers percentage expression ----
  
  #~ pcttable ----
  
  pcttable <- qs::qread(file.path(data_dir, "pcttable.qs"))
  
  stopifnot(identical(
    colnames(pcttable),
    c("gene_name", "id", "pct.exp", "avg.exp", "gene_id")
  ))
  
  
  stopifnot(all(
    pcttable$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    pcttable$gene_id %in% gene_list$gene_id
  ))
  
  
  stopifnot(all(
    pcttable$id %in% all_cell_types
  ))
  
  
  
  
  
  
  
  # Markers precomputed ----
  
  #~ markers ----
  
  markers <- qs::qread(file.path(data_dir, "markers.qs"))
  
  stopifnot(identical(
    colnames(markers) |> sort(),
    c("gene_name", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") |> sort()
  ))
  
  
  stopifnot(all(
    markers$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    markers$gene %in% gene_list$gene_id
  ))
  
  
  stopifnot(all(
    markers$cluster %in% all_neuron_types
  ))
  
  
  
  
  #~ markersAllcells ----
  
  
  markersAllcells <- qs::qread(file.path(data_dir, "markersAllcells.qs"))
  
  stopifnot(identical(
    colnames(markersAllcells) |> sort(),
    c("gene_name", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster") |> sort()
  ))
  
  
  stopifnot(all(
    markersAllcells$gene_name %in% gene_list$gene_name
  ))
  
  stopifnot(all(
    markersAllcells$gene %in% gene_list$gene_id
  ))
  
  
  stopifnot(all(
    markersAllcells$cluster %in% all_cell_types
  ))
  
  
  
  
  
  # Heatmap from file ----
  
  #~ all.TPM.raw ----
  
  all.TPM.raw <- qs::qread(file.path(data_dir, "all.TPM.raw.qs"))
  
  stopifnot(identical(
    colnames(all.TPM.raw),
    all_cell_types
  ))
  
  
  stopifnot(all(
    rownames(all.TPM.raw) %in% gene_list$gene_id
  ))
  
  
  
  
  
  
  #~ TPM.medium ----
  
  TPM.medium <- qs::qread(file.path(data_dir, "TPM.medium.qs"))
  
  stopifnot(identical(
    colnames(TPM.medium),
    all_neuron_types
  ))
  
  
  stopifnot(all(
    rownames(TPM.medium) %in% gene_list$gene_id
  ))
  
  
  
  
  
  
  
  
  
  
  
  
  # DE tests ----
  
  #~ allCells.data ----
  
  allCells.data <- qs::qread(file.path(data_dir, "allCells.data.qs"))
  
  
  stopifnot(all(
    rownames(allCells.data) %in% gene_list$gene_id
  ))
  
  stopifnot(identical(
    unique(TPM.raw.scaled.long$gene_name) |> as.character(),
    gene_list$gene_name[ match(rownames(allCells.data), gene_list$gene_id) ]
  ))
  
  
  #~ allCells.metadata ----
  
  
  allCells.metadata <- qs::qread(file.path(data_dir, "allCells.metadata.qs"))
  
  stopifnot(identical(
    rownames(allCells.metadata), colnames(allCells.data)
  ))
  
  stopifnot(identical(
    unique(allCells.metadata$Cell.type) |> sort(),
    all_cell_types |> sort()
  ))
  
  
  #~ edger_precomputed ----
  
  edger_precomputed <- qs::qread(file.path(data_dir, "edger_precomputed.qs"))
  
  stopifnot(all(
    as.character(unique(edger_precomputed$samples$group)) %in% unique(allCells.metadata$Cell.type)
  ))
  
  
  stopifnot(all(
    rownames(edger_precomputed$counts) %in% gene_list$gene_id
  ))
  
  
  
  
  #~ pseudobulk_matrix ----
  
  pseudobulk_matrix <- qs::qread(file.path(data_dir, "pseudobulk_matrix.qs"))
  
  all.equal(as.matrix(pseudobulk_matrix),
            edger_precomputed$counts)
  
  
  
}



# extra checks for male dataset
validate_male_dataset <- function(dir_male){
  
  # File list ----
  expected_files <- c(
    "comb_allCells.data.qs", "comb_allCells.metadata.qs", "comb_edger_precomputed.qs",
    "comb_pseudobulk_matrix.qs", "herm_all_cell_types.qs"
  )
  
  stopifnot(all(
    expected_files %in% list.files(dir_male)
  ))
  
  all_cell_types_male <- qs::qread(file.path(dir_male, "all_cell_types.qs"))
  all_cell_types_herm <- qs::qread(file.path(dir_male, "herm_all_cell_types.qs"))
  
  gene_list <- qs::qread(file.path(dir_male, "gene_list.qs"))
  
  
  #~ load ----
  
  comb_allCells.data <- qs::qread(file.path(dir_male, "comb_allCells.data.qs"))
  
  comb_allCells.metadata <- qs::qread(file.path(dir_male, "comb_allCells.metadata.qs"))
  
  comb_edger_precomputed <- qs::qread(file.path(dir_male, "comb_edger_precomputed.qs"))
  
  comb_pseudobulk_matrix <- qs::qread(file.path(dir_male, "comb_pseudobulk_matrix.qs"))
  
  
  #~~ genes ----
  
  # comb_allCells.data
  stopifnot(all(
    rownames(comb_allCells.data) %in% gene_list$gene_id
  ))
  
  # comb_pseudobulk_matrix
  stopifnot(all(
    rownames(comb_pseudobulk_matrix) %in% gene_list$gene_id
  ))
  stopifnot(all(
    rownames(comb_pseudobulk_matrix) %in% rownames(comb_allCells.data)
  ))
  
  # comb_edger_precomputed
  stopifnot(all(
    rownames(comb_edger_precomputed$counts) %in% gene_list$gene_id
  ))
  stopifnot(all(
    rownames(comb_edger_precomputed$counts) %in% rownames(comb_allCells.data)
  ))
  
  
  
  #~~ cell names ----
  
  # comb_allCells.data
  stopifnot(identical(
    rownames(comb_allCells.metadata),
    colnames(comb_allCells.data)
  ))
  
  
  # comb_edger_precomputed
  
  # some cell types are filtered out of edgeR object (when too few cells)
  stopifnot(all(
    unique(as.character(comb_edger_precomputed$samples$group))  %in%
          unique(paste0(comb_allCells.metadata$Cell.type, "_", comb_allCells.metadata$Enriched.sex))
  ))
  
  stopifnot(identical(
    colnames(comb_edger_precomputed$counts),
    comb_edger_precomputed$samples$sample_id
  ))
  
  
  # comb_pseudobulk_matrix
  
  stopifnot(all(
    (comb_edger_precomputed$samples$sample_id |>
       str_replace("_Herm_", "_") |>
       str_replace("_Male_", "_")) %in%
      colnames(comb_pseudobulk_matrix)
  ))
  
  
 
  
  stopifnot(all(
    sapply(colnames(comb_pseudobulk_matrix),
           \(cell) any(startsWith(cell, union(all_cell_types_herm, all_cell_types_male))))
  ))
  
  # not true: some cell types have been removed from pseudobulk (namely Unknown_spermatheca_or_uterine_cell, Vulval_muscle_vm2, CP0, HOA)
  # stopifnot(all(
  #   sapply(union(all_cell_types_herm, all_cell_types_male),
  #          \(ct) any(startsWith(colnames(comb_pseudobulk_matrix), ct)))
  # ))
  
  
  
  
}


















