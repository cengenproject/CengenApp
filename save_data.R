
# load("Dataset_1May_2021.rda")

load("Dataset_6July_2021.rda")

# Contains:
#
# allCells     Seurat
# allNeurons   Seurat
# gene_list
# L4.all.TPM.raw           NEW in July 2021
# L4.all.TPM.raw_th        NEW in July 2021
# L4.TPM.medium
# L4.TPM.raw.scaled.long   NEW in July 2021
# markers
# markersAllCells
# med.scaled.long
# pcttable
# ths
# 
# all_cell_types <- sort(unique(allCells$Neuron.type))
# allCells.data <- allCells.data <- GetAssayData(object = allCells[["SCT"]],
#                                                slot = "data")
# allCells.metadata <- allCells@meta.data
# rm(allCells)
# rm(allNeurons)
# 
# save(list = ls(), file = "Dataset_6July_2021_noSeurat.rda")