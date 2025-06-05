
# 2025-03-21 L4 dataset ----


# file.rename("data/2025-03-21/L4.all.TPM.raw_th.qs",
#             "data/2025-03-21/L4.all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th <- qs::qread("data/2025-03-21_L4/L4.all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th |>
#   mutate(threshold = if_else(threshold == "Unfiltered",
#                              "All Cells Unfiltered",
#                              threshold)) |>
#   rename(`Spermathecal-uterine_junction_or_uterine_toroid` = Spermathecal.uterine_junction_or_uterine_toroid,
#          `Unknown_non-neuronal_3` = Unknown_non.neuronal_3) |>
#   qs::qsave("data/2025-03-21_L4/all.TPM.raw_th.qs")
  
  



# file.rename("data/2025-03-21/gene_list.qs",
#             "data/2025-03-21/gene_list.qs.bak")
# 
# gene_list <- qs::qread("data/2025-03-21/gene_list.qs.bak")
# 
# gene_list |>
#   rename(seqnames = seqname) |>
#   qs::qsave("data/2025-03-21/gene_list.qs")



# file.rename("data/2025-03-21_L4/markers.qs",
#             "data/2025-03-21_L4/markers.qs.bak")
# 
# markers <- qs::qread("data/2025-03-21_L4/markers.qs.bak")
# 
# markers |>
#   rename(gene = gene_id) |>
#   qs::qsave("data/2025-03-21_L4/markers.qs")


# file.rename("data/2025-03-21_L4/markersAllcells.qs",
#             "data/2025-03-21_L4/markersAllcells.qs.bak")
# 
# markersAllcells <- qs::qread("data/2025-03-21_L4/markersAllcells.qs.bak")
# 
# markersAllcells |>
#   dplyr::rename(gene = gene_id) |>
#   qs::qsave("data/2025-03-21_L4/markersAllcells.qs")




# file.rename("data/2025-03-21/med.scaled.long.qs",
#             "data/2025-03-21/med.scaled.long.qs.bak")
# 
# med.scaled.long <- qs::qread("data/2025-03-21/med.scaled.long.qs.bak")
# 
# med.scaled.long |>
#   select(-id, -TPM) |>
#   qs::qsave("data/2025-03-21/med.scaled.long.qs")


# file.rename("data/2025-03-21/L4.TPM.raw.scaled.long.qs",
#             "data/2025-03-21/L4.TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long <- qs::qread("data/2025-03-21/L4.TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long |>
#   select(-id, -TPM) |>
#   qs::qsave("data/2025-03-21/TPM.raw.scaled.long.qs")


# file.rename("data/2025-03-21_L4/L4.all.TPM.raw.qs",
#             "data/2025-03-21_L4/all.TPM.raw.qs")
# 
# file.rename("data/2025-03-21_L4/L4.all.TPM.raw_th.qs",
#             "data/2025-03-21_L4/all.TPM.raw_th.qs")
# 
# file.rename("data/2025-03-21_L4/L4.TPM.medium.qs",
#             "data/2025-03-21_L4/TPM.medium.qs")
# 
# file.rename("data/2025-03-21_L4/L4.TPM.raw.scaled.long.qs",
#             "data/2025-03-21_L4/TPM.raw.scaled.long.qs")








# L1 dataset ----



# file.rename("data/2025-03-17_L1/L1.all.TPM.raw.qs",
#             "data/2025-03-17_L1/all.TPM.raw.qs")
# 
# file.rename("data/2025-03-17_L1/L1.all.TPM.raw_th.qs",
#             "data/2025-03-17_L1/all.TPM.raw_th.qs")
# 
# file.rename("data/2025-03-17_L1/L1.TPM.medium.qs",
#             "data/2025-03-17_L1/TPM.medium.qs")
# 
# file.rename("data/2025-03-17_L1/L1.TPM.raw.scaled.long.qs",
#             "data/2025-03-17_L1/TPM.raw.scaled.long.qs")


# file.rename("data/2025-03-17_L1/TPM.raw.scaled.long.qs",
#             "data/2025-03-17_L1/TPM.raw.scaled.long.qs.bak")

# TPM.raw.scaled.long <- qs::qread("data/2025-03-17_L1/TPM.raw.scaled.long.qs.bak")
# TPM.raw.scaled.long |>
#   select(-id, -TPM) |>
#   qs::qsave("data/2025-03-17_L1/TPM.raw.scaled.long.qs")



# file.rename("data/2025-03-17_L1/med.scaled.long.qs",
#             "data/2025-03-17_L1/med.scaled.long.qs.bak")
# 
# med.scaled.long <- qs::qread("data/2025-03-17_L1/med.scaled.long.qs.bak")
# 
# med.scaled.long |>
#   select(-id, -TPM) |>
#   qs::qsave("data/2025-03-17_L1/med.scaled.long.qs")



# file.rename("data/2025-03-17_L1/all.TPM.raw_th.qs",
#             "data/2025-03-17_L1/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th <- qs::qread("data/2025-03-17_L1/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th |>
#   dplyr::rename(
#     `I1?` = `I1.`,
#     `P0.aa/P1.aaa` = `P0.aa.P1.aaa`,
#     `P0.ap/P1.aap` = `P0.ap.P1.aap`,
#     `P1.p?` = `P1.p.`
#   ) |>
#   qs::qsave("data/2025-03-17_L1/all.TPM.raw_th.qs")





# file.rename("data/2025-03-17_L1/ths.qs",
#             "data/2025-03-17_L1/ths.qs.bak")
# 
# ths <- qs::qread("data/2025-03-17_L1/ths.qs.bak")
# 
# ths |>
#   dplyr::rename(
#     `I1?` = `I1.`
#   ) |>
#   qs::qsave("data/2025-03-17_L1/ths.qs")




# file.rename("data/2025-03-17_L1/markers.qs",
#             "data/2025-03-17_L1/markers.qs.bak")
# 
# markers <- qs::qread("data/2025-03-17_L1/markers.qs.bak")
# 
# markers |>
#   dplyr::rename(gene = gene_id) |>
#   qs::qsave("data/2025-03-17_L1/markers.qs")


# file.rename("data/2025-03-17_L1/markersAllcells.qs",
#             "data/2025-03-17_L1/markersAllcells.qs.bak")
# 
# markersAllcells <- qs::qread("data/2025-03-17_L1/markersAllcells.qs.bak")
# 
# markersAllcells |>
#   dplyr::rename(gene = gene_id) |>
#   qs::qsave("data/2025-03-17_L1/markersAllcells.qs")














# adult ----

# ths <- qs::qread("data/010324_adult/ths.qs.bak")
# 
# ths |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/010324_adult/ths.qs")


# file.rename("data/010324_adult/all.TPM.raw_th.qs", "data/010324_adult/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th <- qs::qread("data/010324_adult/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/010324_adult/all.TPM.raw_th.qs")




# file.rename("data/010324_adult/TPM.raw.scaled.long.qs",
#             "data/010324_adult/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long <- qs::qread("data/010324_adult/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long |>
#   dplyr::mutate(tissue = dplyr::if_else(tissue == "Hypodermis",
#                                         "Epidermis", tissue)) |>
#   qs::qsave("data/010324_adult/TPM.raw.scaled.long.qs")






# Male ----

# dir_male <- "data/052225_male/male"
# 
# list.files(dir_male, pattern = "^male_") |>
#   purrr::walk(~ file.rename(
#     file.path(dir_male, .x),
#     file.path(dir_male, stringr::str_remove(.x, "^male_"))
#   ))


# file.rename("data/052225_male/male/all.TPM.raw_th.qs",
#             "data/052225_male/male/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th <- qs::qread("data/052225_male/male/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/052225_male/male/all.TPM.raw_th.qs")



# file.rename("data/052225_male/male/ths.qs",
#             "data/052225_male/male/ths.qs.bak")
# ths <- qs::qread("data/052225_male/male/ths.qs.bak")
# 
# ths |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/052225_male/male/ths.qs")




# file.rename("data/052225_male/male/TPM.raw.scaled.long.qs",
#             "data/052225_male/male/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long <- qs::qread("data/052225_male/male/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long |>
#   dplyr::mutate(tissue = dplyr::if_else(tissue == "Hypodermis",
#                                         "Epidermis", tissue)) |>
#   qs::qsave("data/052225_male/male/TPM.raw.scaled.long.qs")




#~ herm ----

# data_dir <- "data/052225_herm"
# 
# list.files(data_dir, pattern = "^herm") |>
#   purrr::walk(~ file.rename(
#     file.path(data_dir, .x),
#     file.path(data_dir, stringr::str_remove(.x, "^herm[._]"))
#   ))


# file.rename("data/052225_herm/all.TPM.raw_th.qs",
#             "data/052225_herm/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th <- qs::qread("data/052225_herm/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/052225_herm/all.TPM.raw_th.qs")



# file.rename("data/052225_herm/ths.qs",
#             "data/052225_herm/ths.qs.bak")
# ths <- qs::qread("data/052225_herm/ths.qs.bak")
# 
# ths |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/052225_herm/ths.qs")




# file.rename("data/052225_herm/TPM.raw.scaled.long.qs",
#             "data/052225_herm/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long <- qs::qread("data/052225_herm/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long |>
#   dplyr::mutate(tissue = dplyr::if_else(tissue == "Hypodermis",
#                                         "Epidermis", tissue)) |>
#   qs::qsave("data/052225_herm/TPM.raw.scaled.long.qs")




#~ Combine male-herm ----

# male_allCells.metadata <- qs::qread("data/052225_male/male/allCells.metadata.qs")
# herm_allCells.metadata <- qs::qread("data/052225_herm/allCells.metadata.qs")
# 
# male_allCells.metadata$sample_set <- "male"
# herm_allCells.metadata$sample_set <- "herm"
# 
# comb_allCells.metadata <- rbind(
#   male_allCells.metadata,
#   herm_allCells.metadata
# )
# 
# stopifnot(anyDuplicated(rownames(comb_allCells.metadata)) == 0L)
# 
# qs::qsave(comb_allCells.metadata,
#           "data/052225_male/male/comb_allCells.metadata.qs")




# male_allCells.data <- qs::qread("data/052225_male/male/allCells.data.qs")
# herm_allCells.data <- qs::qread("data/052225_herm/allCells.data.qs")
# 
# # match rows (some genes are absent from male, others from herma)
# all_rownames <- union(rownames(male_allCells.data),
#                       rownames(herm_allCells.data))
# 
# 
# rows_missing_from_male <- setdiff(all_rownames, rownames(male_allCells.data))
# rows_missing_from_herm <- setdiff(all_rownames, rownames(herm_allCells.data))
# 
# 
# male_extra_rows <- Matrix::Matrix(0,
#                             nrow = length(rows_missing_from_male),
#                             ncol = ncol(male_allCells.data),
#                             dimnames = list(rows_missing_from_male,
#                                             colnames(male_allCells.data)))
# 
# male_with_zeros <- rbind(male_extra_rows, male_allCells.data)
# male_with_zeros <- male_with_zeros[all_rownames,]
# 
# 
# 
# 
# herm_extra_rows <- Matrix::Matrix(0,
#                                   nrow = length(rows_missing_from_herm),
#                                   ncol = ncol(herm_allCells.data),
#                                   dimnames = list(rows_missing_from_herm,
#                                                   colnames(herm_allCells.data)))
# herm_with_zeros <- rbind(herm_extra_rows, herm_allCells.data)
# herm_with_zeros <- herm_with_zeros[all_rownames,]
# 
# stopifnot(rownames(herm_with_zeros) == rownames(male_with_zeros))
# stopifnot(colnames(herm_with_zeros) == colnames(herm_allCells.data))
# stopifnot(colnames(male_with_zeros) == colnames(male_allCells.data))
# 
# stopifnot(all.equal(
#   colSums(herm_with_zeros),
#   colSums(herm_allCells.data)
# ))
# stopifnot(all.equal(
#   colSums(male_with_zeros),
#   colSums(male_allCells.data)
# ))
# 
# 
# comb_allCells.data <- cbind(
#   male_with_zeros,
#   herm_with_zeros
# )
# 
# stopifnot(anyDuplicated(rownames(comb_allCells.data)) == 0L)
# stopifnot(anyDuplicated(colnames(comb_allCells.data)) == 0L)
# 
# qs::qsave(comb_allCells.data,
#           "data/052225_male/male/comb_allCells.data.qs")


# stopifnot(all.equal(
#   qs::qread("data/052225_male/male/comb_allCells.metadata.qs") |>
#     rownames(),
#   qs::qread("data/052225_male/male/comb_allCells.data.qs") |>
#     colnames()
# ))

                    



                    
# male_pseudobulk_matrix <- qs::qread("data/052225_male/male/pseudobulk_matrix.qs")
# herm_pseudobulk_matrix <- qs::qread("data/052225_herm/pseudobulk_matrix.qs")
# 
# # colnames are different
# stopifnot(length(intersect(
#   colnames(male_pseudobulk_matrix),
#   colnames(herm_pseudobulk_matrix)
# )) == 0L )
# 
# 
# # match rows (some genes are absent from male, others from herma)
# all_rownames <- union(rownames(male_pseudobulk_matrix),
#                       rownames(herm_pseudobulk_matrix))
# 
# 
# rows_missing_from_male <- setdiff(all_rownames, rownames(male_pseudobulk_matrix))
# rows_missing_from_herm <- setdiff(all_rownames, rownames(herm_pseudobulk_matrix))
# 
# 
# male_extra_rows <- Matrix::Matrix(0,
#                                   nrow = length(rows_missing_from_male),
#                                   ncol = ncol(male_pseudobulk_matrix),
#                                   dimnames = list(rows_missing_from_male,
#                                                   colnames(male_pseudobulk_matrix)))
# 
# male_with_zeros <- rbind(male_extra_rows, male_pseudobulk_matrix)
# male_with_zeros <- male_with_zeros[all_rownames,]
# 
# 
# 
# 
# herm_extra_rows <- Matrix::Matrix(0,
#                                   nrow = length(rows_missing_from_herm),
#                                   ncol = ncol(herm_pseudobulk_matrix),
#                                   dimnames = list(rows_missing_from_herm,
#                                                   colnames(herm_pseudobulk_matrix)))
# herm_with_zeros <- rbind(herm_extra_rows, herm_pseudobulk_matrix)
# herm_with_zeros <- herm_with_zeros[all_rownames,]
# 
# stopifnot(rownames(herm_with_zeros) == rownames(male_with_zeros))
# stopifnot(colnames(herm_with_zeros) == colnames(herm_pseudobulk_matrix))
# stopifnot(colnames(male_with_zeros) == colnames(male_pseudobulk_matrix))
# 
# stopifnot(all.equal(
#   colSums(herm_with_zeros),
#   colSums(herm_pseudobulk_matrix)
# ))
# stopifnot(all.equal(
#   colSums(male_with_zeros),
#   colSums(male_pseudobulk_matrix)
# ))
# 
# 
# comb_pseudobulk_matrix <- cbind(
#   male_with_zeros,
#   herm_with_zeros
# )
# 
# attr(comb_pseudobulk_matrix, "sample_set") <- rep(c("male", "herm"),
#                                              times = c(ncol(male_with_zeros),
#                                                        ncol(herm_with_zeros)))
# 
# stopifnot(anyDuplicated(rownames(comb_pseudobulk_matrix)) == 0L)
# stopifnot(anyDuplicated(colnames(comb_pseudobulk_matrix)) == 0L)
# 
# qs::qsave(comb_pseudobulk_matrix,
#           "data/052225_male/male/comb_pseudobulk_matrix.qs")









                    
                    
                    
                    
                    
