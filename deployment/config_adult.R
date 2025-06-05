
# data ----

dataset <- "adult"

data_dir <- "data/052225_herm/"


# aesthetics ----

favicon <- paste0(dataset, ".png")
icon_big <- paste0(dataset, "_big.png")
bg_color <-"#688970"



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

paper_authors <- "Claire Olson, Seth Taylor, Rebecca McWhirter, Brooks Leyhew, Chen Wang, Giulio Valperga, Isabel Courtney, Marc Hammarlund, Oliver Hobert, David M. Miller, III."
paper_title <- "Sex-specific gene expression for an adult nervous system,"
paper_date_venue <- "manuscript in preparation."
paper_footnote <- NULL


