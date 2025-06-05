
# data ----

dataset <- "male"

data_dir <- "data/052225_male/male/"

herm_all_cell_types_path <- "data/052225_male/male/herm_all_cell_types.qs"



# aesthetics ----

favicon <- paste0(dataset, ".png")
icon_big <- paste0(dataset, "_big.png")
bg_color <- "#9f8200"



# Warning message ----
# May contain a list of unreliable genes and/or a warning about methanol

# library(wbData)
# gids <- wb_load_gene_ids(295)
# unreliable_gene_ids <- c()
# unreliable_gene_ids |> i2s(gids, warn_missing = TRUE) |> paste(collapse = "', '")
unreliable_gene_names <- c()


info_msg <- div(
  "Note: this data was generated from methanol-fixed cells. Expression profiles might differ from previous datasets (generated from live cells)."
  ,
  class = "alert alert-warning alert-light"
)


# paper reference ----

paper_authors <- "Olson, et al."
paper_title <- "Sex-specific gene expression for an adult nervous system,"
paper_date_venue <- "manuscript in preparation."
paper_footnote <- NULL



