
# data ----

dataset <- "L1"

data_dir <- "data/2025-03-17_L1"


# aesthetics ----

favicon <- "L1.png"
icon_big <- "L1_big.png"
bg_color <- "#b4d0f9"



# Warning message ----
# May contain a list of unreliable genes and/or a warning about methanol

# library(wbData)
# gids <- wb_load_gene_ids(295)
# unreliable_gene_ids <- c(
#   "WBGene00015522", "WBGene00196150", "WBGene00206533", "WBGene00000429", "WBGene00001079", "WBGene00001135", "WBGene00018172",
#   "WBGene00200518", "WBGene00015520", "WBGene00015521", "WBGene00023498", "WBGene00023497", "WBGene00004397", "WBGene00006843",
#   "WBGene00006783", "WBGene00006789"
# )
# unreliable_gene_ids |> i2s(gids, warn_missing = TRUE) |> paste(collapse = "', '")
unreliable_gene_names <- c('C06E1.7', 'C30A5.11', 'C30A5.16', 'ceh-2', 'dpy-20',
                           'eat-4', 'F38B6.2', 'F38B6.17', 'fip-3', 'fipr-16', 'lin-15A',
                           'lin-15B', 'rol-6', 'unc-119', 'unc-47', 'unc-54')


info_msg <- div(
    "WARNING: Expression values for ",
    unreliable_gene_names |> paste(collapse = ", "),
    " are unreliable as they have been overexpressed to generate transgenic strains."
  ,
  class = "alert alert-warning alert-light"
)


# paper reference ----

paper_authors <- "Seth R. Taylor*, Claire Olson, Rebecca McWhirter, Sidharth Goel, Lidia Ripoll-Sanchez, S Talmage Barney, Drew Hardin, Alexis Rolfson, Jacob Pattee, Alexander Atkinson, Ethan Grundvig, Isabel Courtney, Giulio Valperga, G. Robert Aguilar, Daniel M. Merrit, Isabel Beets, Petra Vertes, William R. Schafer, Erdem Varol, Marc Hammarlund, Oliver Hobert, David M. Miller, III"
paper_title <- "A gene expression atlas of neurogenic cells in the C. elegans first stage larva,"
paper_date_venue <- "manuscript in preparation."
paper_footnote <- "*Corresponding author."


