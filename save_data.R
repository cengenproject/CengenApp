# # Load the full datasets in the Shiny app, and save it without the Seurat objects.
# Compute some other objects which might be useful.


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


# For sc Wilcoxon tests
all_cell_types <- sort(unique(allCells$Neuron.type))
allCells.data <- allCells.data <- GetAssayData(object = allCells[["SCT"]],
                                               slot = "data")
allCells.metadata <- allCells@meta.data




# For pseudobulk tests
pseudobulk_matrix <- AggregateExpression(allCells,
                                  assays = "RNA",
                                  slot = "counts",
                                  group.by = c("Neuron.type", "SampleID"))[["RNA"]]

pseudosamples <- as.character(colnames(pseudobulk_matrix))
regx <- stringr::str_match(pseudosamples,
                           "^([A-Za-z0-9\\-_]+)_([0-9]{4}-ST-[12])$")
pseudobulk_metadata <- data.frame(sample = pseudosamples,
                                  cell_type = regx[,2],
                                  batch = regx[,3])

rm(regx); rm(pseudosamples)

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

to_save <- c("gene_list",
             "L4.all.TPM.raw",
             "L4.all.TPM.raw_th",
             "L4.TPM.medium",
             "L4.TPM.raw.scaled.long",
             "markers",
             "markersAllcells",
             "med.scaled.long",
             "pcttable",
             "ths",
             "all_cell_types",
             "allCells.data",
             "allCells.metadata",
             "pseudobulk_matrix",
             "pseudobulk_metadata",
             "edger_precomputed")



save(list = to_save, file = "Dataset_6July_2021_noSeurat2.rda")




