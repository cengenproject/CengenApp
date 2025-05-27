
# 2025-03-21 L4 dataset ----


# file.rename("data/2025-03-21/L4.all.TPM.raw_th.qs",
#             "data/2025-03-21/L4.all.TPM.raw_th.qs.bak")
# 
# L4.all.TPM.raw_th <- qs::qread("data/2025-03-21/L4.all.TPM.raw_th.qs.bak")
# 
# L4.all.TPM.raw_th |>
#   mutate(threshold = if_else(threshold == "Unfiltered",
#                              "All Cells Unfiltered",
#                              threshold)) |>
#   qs::qsave("data/2025-03-21/L4.all.TPM.raw_th.qs")


# file.rename("data/2025-03-21/gene_list.qs",
#             "data/2025-03-21/gene_list.qs.bak")
# 
# gene_list <- qs::qread("data/2025-03-21/gene_list.qs.bak")
# 
# gene_list |>
#   rename(seqnames = seqname) |>
#   qs::qsave("data/2025-03-21/gene_list.qs")



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


# Male ----

# dir_male <- "data/052225_male/male"
# 
# list.files(dir_male, pattern = "^male_") |>
#   purrr::walk(~ file.rename(
#     file.path(dir_male, .x),
#     file.path(dir_male, stringr::str_remove(.x, "^male_"))
#   ))






