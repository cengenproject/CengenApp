
# read config file
config <- list.files("deployment/", pattern = "config_", full.names = TRUE)

if(length(config) == 1){
  source(config)
} else if(length(config) > 1){
  
  # here set default file
  source("deployment/config_male.R")
  
} else{
  stop("Configuration file not found")
}




message("-------------  dataset: ", dataset," ----------------")

# check config loaded
stopifnot(nchar(dataset) > 1)
stopifnot(file.exists(file.path("www", favicon) ))
stopifnot(file.exists(file.path("www", icon_big) ))
stopifnot(dir.exists(data_dir))
stopifnot( length(list.files(data_dir)) > 5 )



# variables directly defined from config file

favicon <- paste0(dataset, ".png")
icon_big <- paste0(dataset, "_big.png")


available_apps <- data.frame(
  name = c("L1", "L4", "adult", "male"),
  url = c("L1app", "L4app", "adult", "male")
)

other_apps <- available_apps[available_apps$name != dataset,]




# global variables for all apps

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



