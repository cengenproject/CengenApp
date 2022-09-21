# library(readr)
# library(dplyr)
# library(wbData)
# 
# 
# 
# 
# load("Dataset_1May_2021.rda")
# 
# # Contains:
# # 
# # allCells     Seurat
# # allNeurons   Seurat
# # gene_list
# # L4.TPM.medium
# # markers
# # markersAllCells
# # med.scaled.long
# # pcttable
# # ths
# 
# all_cell_types <- sort(unique(allCells$Neuron.type))
# allCells.data <- allCells.data <- GetAssayData(object = allCells[["SCT"]], 
#                                                slot = "data")
# allCells.metadata <- allCells@meta.data
# rm(allCells)
# rm(allNeurons)
# 
# save(list = ls(), file = "Dataset_1May_2021_noSeurat.rda")
# 
# load("Dataset_1May_2021_noSeurat.rda")
# 
# # save as DuckDB
# con <- DBI::dbConnect(duckdb::duckdb(), "data/t_exp.duckdb.db")
# DBI::dbWriteTable(con, name = "t_exp",value = tx_long)
# 
# 
# 
# # to read
# t_exp_db <- pool::dbPool(RSQLite::SQLite(),
#                          dbname = "data/t_exp.sqlite.db",
#                          read_only = TRUE)
# 
# onStop(function() {
#   pool::poolClose(t_exp_db)
# })