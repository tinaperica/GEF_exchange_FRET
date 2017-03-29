### GTPase assay dirty water script

reaction.V <- 80 # ul, for reactions in 96-well plate
buffer_with_Mg.V <- reaction.V / 2
get_volume <- function (stock, final) {
  round( reaction.V / (stock/final) , 1 ) 
}
proteins <- c("PE1:GDP")
#### stock concentrations
stocks <- list()
stocks[["PE1:GDP"]] <- c(12.67) # uM
stocks[["PE1:GTP"]] <- c(8.84) # uM
stocks[["sensor"]] <- 10  # original stock is 500 uM, diluted in the assay buffer 10x
stocks[["GAP"]] <- 46.95 # uM
stocks[["GEF"]] <- 101.32/2 # uM
stocks[["MOG1"]] <- 596/50 # uM
stocks[["GTP"]] <- 300 # uM
### final reaction concentrations
final_c <- list()
final_c[["PE1:GDP"]] <- c(2)
#final_c[["PE1:GTP"]] <- c(2) # uM
final_c[["sensor"]] <- c(2) # uM
final_c[["GAP"]] <- c(0, 0.01) # uM
final_c[["GEF"]] <- c(0, 0.005) # uM
#final_c[["GTP"]] <- c() # uM
final_c[["MOG1"]] <- c(0, 0.01, 0.1, 2)
experiment_combinations <- expand.grid( final_c )
final <- data.frame()
for (protein in proteins) {
  experiment_volumes <- list()
  for ( i in seq_along(names(final_c)) ) {
    component <- names(final_c)[i]
    volume <- get_volume( stocks[[component]], final_c[[component]] )
    if ( length(volume) > 0) {
      experiment_volumes[[paste0(component, ".volume")]] <- volume
    }
  }
  exp_volume_combinations <- expand.grid(experiment_volumes)
  exp_volume_combinations <- exp_volume_combinations[, !names(exp_volume_combinations) %in% c("GAP.volume", "GEF.volume")]
  exp_volume_combinations <- exp_volume_combinations[ rowSums(exp_volume_combinations == 0) != ncol(exp_volume_combinations), ]  ### this removes the rows with all 0.0
  buffer_wo_Mg.V <- (reaction.V/2) - rowSums( exp_volume_combinations )
  exp_comb <- cbind (cbind(exp_volume_combinations, experiment_combinations), data.frame( "buffer_wo_Mg_V" = buffer_wo_Mg.V, "enzyme_mix_with_Mg" = buffer_with_Mg.V) )
  final <- rbind( final, cbind(exp_comb, "protein" = protein) )
}

colSums(final[, 1:5])
final

