### dirty water script for GTP loading

### GTPase assay dirty water script
reaction_volume <- 250 # ul
buffer_components <- list("name" = c("EDTA", "Tris", "NaCl", "DTT"),
                          "stock" = c(100, 1000, 3000, 100), 
                          "final_conc" = c(10, 50, 100, 1)
                      )
buffer_components[["final_volumes"]] <- reaction_volume / (buffer_components[["stock"]]/buffer_components[["final_conc"]])

GTP_stock <- 90 # mM conc in nmol/ul
desired_GTP_to_protein_ratio <- 200
protein_stocks <- read.delim("current_protein_stocks_in_mg_per_ml.txt", head = F)
names(protein_stocks) <- expression(number, mutant, concentration)
protein_stocks[["protein"]] <- paste0("PE", protein_stocks$number)
Gsp1_MW <- 24.8 # kg/mol
protein_stocks[["mM_conc"]] <- protein_stocks$concentration/Gsp1_MW  # nmol/ul
volume_for_protein_and_GTP <- reaction_volume - sum(buffer_components[["final_volumes"]])
volume_for_protein_and_GTP # 171 ul

get_protein_nmol <- function(protein_stock_mM) {
  protein_nmol <- (volume_for_protein_and_GTP * GTP_stock * protein_stock_mM) / (protein_stock_mM * desired_GTP_to_protein_ratio + GTP_stock) 
  return(protein_nmol)
}
loading_table <- data.frame()
for ( i in seq_along(protein_stocks[["protein"]]) ) {
  protein <- protein_stocks[["protein"]][i]
  protein_stock_mM <- protein_stocks$mM_conc[ protein_stocks[["protein"]] == protein]
  protein_nmol <- round( get_protein_nmol(protein_stock_mM), 2)
  protein_V <- round(protein_nmol / protein_stock_mM, 2)
  GTP_V <- round(volume_for_protein_and_GTP - protein_V, 2)
  GTP_nmol <- round( GTP_stock * GTP_V, 2)
  GTP_protein_ratio <- round(GTP_nmol / protein_nmol,0)
  loading_table <- rbind(loading_table, data.frame(protein, protein_nmol, protein_V, 
                                                    GTP_nmol, GTP_V,
                                                   "EDTA" = buffer_components[["final_volumes"]][1],
                                                   "Tris" = buffer_components[["final_volumes"]][2],
                                                   "NaCl" = buffer_components[["final_volumes"]][3],
                                                   "DTT" = buffer_components[["final_volumes"]][4],
                                                   GTP_protein_ratio))
}
loading_table
