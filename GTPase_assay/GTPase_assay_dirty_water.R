### GTPase assay dirty water script

reaction.volume <- 100 # ul, for reactions in 96-well plate
buffer_with_Mg.V <- reaction.volume / 2
get_volume <- function (stock, final) {
  round( reaction.volume / (stock/final) , 1 ) 
}
proteins <- c( 
              "PE1:GTP_1", "PE1:GTP_2"
              )
protein_stocks <- c(
  25.41, 126.05
              ) # uM
sensor_stock <- 500 # original stock is 500 uM, diluted in the assay buffer 5x
#GAP_stock <- 5 # uM
### final reaction concentrations of Gsp1:GTP (order the same as in the proteins vector)
final_c <- list()
final_c[[1]] <- c() # uM
final_c[[2]] <- c(4) # uM

#### Gsp1 concentrations to run without GAP (to measure intrinsinc hydrolysis)
intrinsinc_c <- list()
intrinsinc_c[[1]] <- c(1)

final_sensor <- 18 # uM
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


