
# note the assignment in the global environment
load_as_needed <- function(dataset){
  if(!exists(dataset)){
    assign(dataset,
           qs::qread(file.path(data_dir, paste0(dataset, ".qs"))),
           envir = .GlobalEnv)
  }
}
