unreliable_gene_ids <- c(
  "WBGene00023498",
  "WBGene00023497",
  "WBGene00004397",
  "WBGene00006843",
  "WBGene00004010",
  "WBGene00006789",
  "WBGene00001135",
  "WBGene00001079",
  "WBGene00006783",
  "WBGene00000501",
  "WBGene00006788",
  "WBGene00001555",
  "WBGene00206533",
  "WBGene00011964",
  "WBGene00018172",
  "WBGene00016259",
  "WBGene00023407"
)

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
