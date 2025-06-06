
# if there are several configuration files available, which one to use as default
default_config <- "deployment/config_adult.R"


# shared variables for all apps

ordered_tissues <- c(
  "Neuron", "Neuronal_progenitors", "Glia", "Epidermis", "Muscle_mesoderm", "Pharynx",
  "Excretory","Intestine","Rectal_cells", "Reproductive",
  "Unknown", "Unannotated"
)

pharyngeal_neurons <- c(
  paste0("I", 1:6),
  paste0("M", 1:5 |> c("C","I")),
  "NSM"
)

male_neurons <- c("CA1-4", "CA5-6", "CA7", "CA8-9", "CEM", "CP0", "CP1-4", "CP5-6",
                  "CP7-8", "CP9", "DVE_DVF", "DX", "EF", "HOA", "HOB?", "male.44A",
                  "male.44B", "male.45", "male.46", "male.47", "male.48", "male.49",
                  "male.53", "male.55A", "male.55B", "male.56", "MCM", "PCA?", "PCC",
                  "PDC", "PGA", "PHD?", "PVV", "PVX", "PVY", "R1B_R4B_R7B?", "R2B?",
                  "R3B_R9B?", "R4A", "R5A", "R5B?", "R6B", "R7A", "R8A?", "R8B?",
                  "R9A", "SPD", "SPV")




#--- no need to update below ---


# read config file
config <- list.files("deployment/", pattern = "config_", full.names = TRUE)

if(length(config) == 1){
  
  source(config)
} else if(length(config) > 1){
  
  source(default_config)
} else{
  
  stop("Configuration file not found")
}
rm(config)



message("-------------  dataset: ", dataset," ----------------")

# check config loaded
stopifnot(nchar(dataset) > 1)
stopifnot(dir.exists(data_dir))
stopifnot( length(list.files(data_dir)) > 5 )
stopifnot(nchar(bg_color) > 4)

# load global data and functions

source("R/load_as_needed.R")
source("R/perform_de.R")
source("R/plot_heatmap.R")

load_as_needed("gene_list")
load_as_needed("all_cell_types")
load_as_needed("all_neuron_types")

if(dataset == "male"){
  male_all_cell_types <- qs::qread(file.path(data_dir, "all_cell_types.qs"))
  herm_all_cell_types <- qs::qread(herm_all_cell_types_path)
}



# variables directly defined from config file

favicon <- paste0(dataset, ".png")
icon_big <- paste0(dataset, "_big.png")


available_apps <- data.frame(
  name = c("L1", "L4", "adult", "male"),
  url = c("L1app", "L4app", "adult", "male")
)

other_apps <- available_apps[available_apps$name != dataset,]


footer <- div(
  h5("Citation:"),
  paper_authors, br(),
  em(paper_title),
  paper_date_venue, br(),
  paper_footnote,
  class = "alert alert-info mt-5"
)





