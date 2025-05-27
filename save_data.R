
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




# herm ----

# data_dir <- "data/052225_male/herm/"
# 
# list.files(data_dir, pattern = "^herm") |>
#   purrr::walk(~ file.rename(
#     file.path(data_dir, .x),
#     file.path(data_dir, stringr::str_remove(.x, "^herm[._]"))
#   ))


# file.rename("data/052225_male/herm/all.TPM.raw_th.qs",
#             "data/052225_male/herm/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th <- qs::qread("data/052225_male/herm/all.TPM.raw_th.qs.bak")
# 
# all.TPM.raw_th |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/052225_male/herm/all.TPM.raw_th.qs")



# file.rename("data/052225_male/herm/ths.qs",
#             "data/052225_male/herm/ths.qs.bak")
# ths <- qs::qread("data/052225_male/herm/ths.qs.bak")
# 
# ths |>
#   dplyr::rename(id = X) |>
#   qs::qsave("data/052225_male/herm/ths.qs")




# file.rename("data/052225_male/herm/TPM.raw.scaled.long.qs",
#             "data/052225_male/herm/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long <- qs::qread("data/052225_male/herm/TPM.raw.scaled.long.qs.bak")
# 
# TPM.raw.scaled.long |>
#   dplyr::mutate(tissue = dplyr::if_else(tissue == "Hypodermis",
#                                         "Epidermis", tissue)) |>
#   qs::qsave("data/052225_male/herm/TPM.raw.scaled.long.qs")




