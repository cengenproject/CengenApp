
# data ----
dataset <- "L4"

data_dir <- "data/2025-03-21_L4"

# aesthetics ----
bg_color <- "#bc9bac"



# Warning message ----
# May contain a list of unreliable genes and/or a warning about methanol

# library(wbData)
# gids <- wb_load_gene_ids(295)
# unreliable_gene_ids <- c(
#   "WBGene00023498",
#   "WBGene00023497",
#   "WBGene00004397",
#   "WBGene00006843",
#   "WBGene00004010",
#   "WBGene00006789",
#   "WBGene00001135",
#   "WBGene00001079",
#   "WBGene00006783",
#   "WBGene00000501",
#   "WBGene00006788",
#   "WBGene00001555",
#   "WBGene00206533",
#   "WBGene00011964",
#   "WBGene00018172",
#   "WBGene00016259",
#   "WBGene00023407"
# )
# unreliable_gene_ids |> i2s(gids, warn_missing = TRUE) |> paste(collapse = "', '")
unreliable_gene_names <- c('lin-15A', 'lin-15B', 'rol-6', 'unc-119', 'pha-1', 'unc-54',
                           'eat-4', 'dpy-20', 'unc-47', 'cho-1', 'unc-53', 'gcy-35',
                           'C30A5.16', 'saeg-2', 'F38B6.2', 'C30F8.3', 'cex-1')


info_msg <- div(
  "WARNING: Expression values for ",
  unreliable_gene_names |> paste(collapse = ", "),
  " are unreliable as they have been overexpressed to generate transgenic strains."
  ,
  class = "alert alert-warning alert-light"
)



# paper reference ----

paper_authors <- "Seth R. Taylor, Gabriel Santpere, Alexis Weinreb, Alec Barrett, Molly B. Reilly, Chuan Xu, Erdem Varol, Panos Oikonomou, Lori Glenwinkel, Rebecca McWhirter, Abigail Poff, Manasa Basavaraju, Ibnul Rafi, Eviatar Yemini, Steven J. Cook, Alexander Abrams, Berta Vidal, Cyril Cros, Saeed Tavazoie, Nenad Sestan, Marc Hammarlund, Oliver Hobert, David M. Miller III."
paper_title <- "Molecular topography of an entire nervous system,"
paper_date_venue <- "Cell (2021)."
paper_footnote <- a("https://doi.org/10.1016/j.cell.2021.06.023", href = "https://doi.org/10.1016/j.cell.2021.06.023")
