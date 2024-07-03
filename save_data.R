# # Load the full datasets in the Shiny app, and save it without the Seurat objects.
# Compute some other objects which might be useful.

library(Seurat)
load("Dataset_6July_2021.rda")

# Contains:
#
# allCells     Seurat -> to remove
# allNeurons   Seurat -> to remove
# gene_list
# L4.all.TPM.raw           NEW in July 2021
# L4.all.TPM.raw_th        NEW in July 2021
# L4.TPM.medium
# L4.TPM.raw.scaled.long   NEW in July 2021
# markers
# markersAllcells
# med.scaled.long
# pcttable
# ths


# general
all_cell_types <- sort(unique(allCells$Neuron.type))
all_neuron_types <- colnames(L4.TPM.medium)


# correct genelist (WBGene00007396 is as tiar-3 in `med.scaled.long` and `gene.list` and as rnp-9 in `L4.TPM.raw.scaled.long`)
L4.TPM.raw.scaled.long <- qs::qread("data/L4.TPM.raw.scaled.long.qs.bak")
levels(L4.TPM.raw.scaled.long$gene_name)[levels(L4.TPM.raw.scaled.long$gene_name) == "rnp-9"] <- "tiar-3"
levels(L4.TPM.raw.scaled.long$gene_name) <- gsub("rnp-9","tiar-3", levels(L4.TPM.raw.scaled.long$gene_name))

# qs::qsave(L4.TPM.raw.scaled.long, "data/L4.TPM.raw.scaled.long.qs")


# For sc Wilcoxon tests
allCells.data <- allCells.data <- GetAssayData(object = allCells[["SCT"]],
                                               slot = "data")
allCells.metadata <- allCells@meta.data



# For marker tables

markers$p_val <-
  formatC(markers$p_val, format = "e", digits = 3) %>% gsub(" ", "", .)
markers$p_val_adj <-
  formatC(markers$p_val_adj, format = "e", digits = 3) %>% gsub(" ", "", .)
markers$avg_log2FC <-
  formatC(markers$avg_log2FC, digits = 3) %>% gsub(" ", "", .)
markersAllcells$p_val <-
  formatC(markersAllcells$p_val,
          format = "e",
          digits = 3) %>% gsub(" ", "", .)
markersAllcells$p_val_adj <-
  formatC(markersAllcells$p_val_adj,
          format = "e",
          digits = 3) %>% gsub(" ", "", .)
markersAllcells$avg_log2FC <-
  formatC(markersAllcells$avg_log2FC, digits = 3) %>% gsub(" ", "", .)




# For pseudobulk tests
pseudobulk_matrix <- AggregateExpression(allCells,
                                  assays = "RNA",
                                  slot = "counts",
                                  group.by = c("Neuron.type", "SampleID"))[["RNA"]]

pseudosample_counts <- allCells@meta.data |>
  count(Neuron.type, SampleID) |>
  mutate(sample_id = paste(Neuron.type, SampleID, sep = "_"))


stopifnot(all.equal(pseudosample_counts$sample_id,
                    as.character(colnames(pseudobulk_matrix))))

# filter out replicates with <10 single cells
pseudobulk_matrix <- pseudobulk_matrix[,pseudosample_counts$sample_id[pseudosample_counts$n > 10]]

pseudobulk_metadata <- pseudosample_counts |>
  filter(n > 10) |>
  rename(sample_id = sample_id,
         cell_type = Neuron.type,
         batch = SampleID,
         nb_single_cells = n)


# Precompute edgeR object
library(edgeR)
edger_precomputed <- DGEList(counts=pseudobulk_matrix,
             samples = pseudobulk_metadata,
             group = pseudobulk_metadata$cell_type)


keep <- filterByExpr(edger_precomputed)
edger_precomputed <- edger_precomputed[keep, , keep.lib.sizes=FALSE]
edger_precomputed <- calcNormFactors(edger_precomputed)

edger_precomputed <- estimateDisp(edger_precomputed)

# overwrite to use for pseudobulk Wilcoxon test
pseudobulk_matrix <- pseudobulk_matrix[keep,]



rm(allCells)
rm(allNeurons)

# to_save <- c("gene_list",
#              "L4.all.TPM.raw",
#              "L4.all.TPM.raw_th",
#              "L4.TPM.medium",
#              "L4.TPM.raw.scaled.long",
#              "markers",
#              "markersAllcells",
#              "med.scaled.long",
#              "pcttable",
#              "ths",
#              "all_cell_types",
#              "allCells.data",
#              "allCells.metadata",
#              "pseudobulk_matrix",
#              "pseudobulk_metadata",
#              "edger_precomputed")
# 
# 
# 
# save(list = to_save, file = "Dataset_6July_2021_noSeurat2.rda")


## Save data in separate files to load only those needed
library(qs)

qsave(gene_list, "data/gene_list.qs")
qsave(L4.all.TPM.raw, "data/L4.all.TPM.raw.qs")
qsave(L4.all.TPM.raw_th, "data/L4.all.TPM.raw_th.qs")
qsave(L4.TPM.medium, "data/L4.TPM.medium.qs")
qsave(L4.TPM.raw.scaled.long, "data/L4.TPM.raw.scaled.long.qs")
qsave(markers, "data/markers.qs")
qsave(markersAllcells, "data/markersAllcells.qs")
qsave(med.scaled.long, "data/med.scaled.long.qs")
qsave(pcttable, "data/pcttable.qs")
qsave(ths, "data/ths.qs")
qsave(all_cell_types, "data/all_cell_types.qs")
qsave(allCells.data, "data/allCells.data.qs")
qsave(allCells.metadata, "data/allCells.metadata.qs")
qsave(pseudobulk_matrix, "data/pseudobulk_matrix.qs")
# qsave(pseudobulk_metadata, "data/pseudobulk_metadata.qs")
qsave(edger_precomputed, "data/edger_precomputed.qs")
qsave(all_neuron_types, "data/all_neuron_types.qs")

