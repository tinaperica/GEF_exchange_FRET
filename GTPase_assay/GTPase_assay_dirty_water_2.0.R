### GTPase assay dirty water script

reaction.volume <- 100 # ul, for reactions in 96-well plate
buffer_with_Mg.V <- reaction.volume / 2
get_volume <- function (stock, final) {
  round( reaction.volume / (stock/final) , 1 ) 
}
proteins <- c( 
  "PE1:GTP_5x_diluted", "PE1:GTP")
protein_stocks <- c(
          14.17, 70.85) # uM
sensor_stock <- vector()
sensor_stock[1] <- 100 # original stock is 500 uM, diluted in the assay buffer 5x
#sensor_stock[2] <- 250 ### original stock is 500 uM, diluted in the assay buffer 2x
### final reaction concentrations of Gsp1:GTP (order the same as in the proteins vector)
final_c <- list()
final_c[[1]] <- c(0.5, 1)
final_c[[2]] <- c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10)

#### Gsp1 concentrations to run without GAP (to measure intrinsinc hydrolysis)
intrinsinc_c <- list()
intrinsic_c[[1]] <- c(3)

#final_sensor <- 18 # uM
#sensor.V = get_volume(sensor_stock, final_sensor)
final_GAP <- 15 # nM

reaction.mix.table <- data.frame()
for (i in seq_along(proteins)) {
  prot_sample <- proteins[i]
  prot <- strsplit(prot_sample, ":")
  stock.conc <- protein_stocks[i]
  for (j in seq_along( final_c[[i]]) ) {
    Ran.conc <- final_c[[i]][j]
    if (Ran.conc <= 2) {
      final_sensor <- 4
      gain <- 74
      sensor_stock_to_use <- sensor_stock[1]
      sensor.V <- get_volume(sensor_stock_to_use, final_sensor)
    } else if (Ran.conc > 2 & Ran.conc <= 6) {
      final_sensor <- 10
      gain <- 74
      sensor_stock_to_use <- sensor_stock[1]
      sensor.V <- get_volume(sensor_stock_to_use, final_sensor)
    } else if (Ran.conc > 6) {
      final_sensor <- 18
      gain <- 74
      sensor_stock_to_use <- sensor_stock[1]
      sensor.V <- get_volume(sensor_stock_to_use, final_sensor)
    }
    Ran.volume <- get_volume(stock.conc, Ran.conc)
    buffer_wo_Mg.volume = round((reaction.volume/2) - (Ran.volume + sensor.V), 2)
    reaction.mix.table <- rbind(reaction.mix.table,
        data.frame("protein" = prot[[1]][1], prot_sample, Ran.conc, Ran.volume, buffer_wo_Mg.volume,
                   final_sensor, sensor_stock_to_use, sensor.V, "GAP" = final_GAP, gain))
    
  }
  if (length(intrinsinc_c) >= i) {
    for (k in seq_along( intrinsinc_c[[i]] )) {
      Ran.conc <- intrinsinc_c[[i]][k]
      Ran.volume <- get_volume(stock.conc, Ran.conc)
      buffer_wo_Mg.volume = round((reaction.volume/2) - (Ran.volume + sensor.V), 2)
      reaction.mix.table <- rbind(reaction.mix.table,
          data.frame("protein" = prot[[1]][1], prot_sample, Ran.conc, Ran.volume, buffer_wo_Mg.volume,
                     final_sensor, sensor_stock_to_use, sensor.V, "GAP" = 0, gain))
    }
  }
}
reaction.mix.table <- reaction.mix.table[order(
  reaction.mix.table$final_sensor, reaction.mix.table$protein, reaction.mix.table$Ran.conc), ]
reaction.mix.table

sum(reaction.mix.table$buffer_wo_Mg.volume)
sum(reaction.mix.table$sensor.V)
with(reaction.mix.table, aggregate(Ran.volume, by = list(prot_sample = prot_sample), sum))


