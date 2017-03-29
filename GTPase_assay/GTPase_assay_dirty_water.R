### GTPase assay dirty water script

reaction.volume <- 100 # ul, for reactions in 96-well plate
buffer_with_Mg.V <- reaction.volume / 2
get_volume <- function (stock, final) {
  round( reaction.volume / (stock/final) , 1 ) 
}
proteins <- c( 
              "PE3:GTP_LOW", "PE3:GTP_HIGH",
              "PE4:GTP_LOW", "PE4:GTP_HIGH",
              
              "PE13:GTP_LOW", "PE13:GTP_HIGH"
              )
protein_stocks <- c(
                    10.44, 17.42,
                    19.49, 29.55,
                    
                    28.60, 40.53
                    ) # uM
sensor_stock <- 100  # original stock is 500 uM, diluted in the assay buffer 5x
#GAP_stock <- 5 # uM
### final reaction concentrations of Gsp1:GTP (order the same as in the proteins vector)
final_c <- list()
#final_c[[1]] <- c(2) # uM
#final_c[[2]] <- c(3) # uM
final_c[[3]] <- c(1, 3, 4) # uM
final_c[[4]] <- c(2) # uM
final_c[[5]] <- c(1, 2, 3) # uM
final_c[[6]] <- c(3, 4, 5, 6, 7) # uM


#### Gsp1 concentrations to run without GAP (to measure intrinsinc hydrolysis)
intrinsinc_c <- list()
# intrinsinc_c[[2]] <- c(3)
# intrinsinc_c[[3]] <- c(3)
# intrinsinc_c[[4]] <- c(3)
# intrinsinc_c[[6]] <- c(3)
# intrinsinc_c[[7]] <- c(3)
# intrinsinc_c[[8]] <- c(3)

final_sensor <- 4 # uM
sensor.V = get_volume(sensor_stock, final_sensor)
final_GAP <- 6 # nM

reaction.mix.table <- data.frame()
for (i in seq_along(proteins)) {
  prot <- proteins[i]
  stock.conc <- protein_stocks[i]
  for (j in seq_along( final_c[[i]]) ) {
    Ran.conc <- final_c[[i]][j]
    Ran.volume <- get_volume(stock.conc, Ran.conc)
    buffer_wo_Mg.volume = round((reaction.volume/2) - (Ran.volume + sensor.V), 2)
    reaction.mix.table <- rbind(reaction.mix.table,
        data.frame(prot, Ran.conc, Ran.volume, buffer_wo_Mg.volume, sensor.V, "GAP" = final_GAP))
    
  }
  if (length(intrinsinc_c) >= i) {
    for (k in seq_along( intrinsinc_c[[i]] )) {
      Ran.conc <- intrinsinc_c[[i]][k]
      Ran.volume <- get_volume(stock.conc, Ran.conc)
      buffer_wo_Mg.volume = round((reaction.volume/2) - (Ran.volume + sensor.V), 2)
      reaction.mix.table <- rbind(reaction.mix.table,
          data.frame(prot, Ran.conc, Ran.volume, buffer_wo_Mg.volume, sensor.V, "GAP" = 0))
    }
  }
}

reaction.mix.table

sum(reaction.mix.table$buffer_wo_Mg.volume)
sum(reaction.mix.table$sensor.V)
with(reaction.mix.table, aggregate(Ran.volume, by = list(prot = prot), sum))


