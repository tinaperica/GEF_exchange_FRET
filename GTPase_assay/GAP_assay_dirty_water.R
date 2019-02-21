library(tidyverse)
SubstrateDilute <- function(stock, work, volume, Mw) {
  molar.stock <- 1000*(stock/Mw)
  dilutionX <- molar.stock/work
  protein.V <- round(volume/dilutionX,1)
  buffer.V <- round(volume - protein.V,1)
  return(data.frame(protein.V, buffer.V, "dilution" = round(dilutionX, 2)))
}

### Set up  experiment parameters
plate <- 7
row <- "E"
date <- "20190220"
Rans <- c("PE3_T34E")
#sensor.working <- 10 # uM

outdir <- "GTPase_assay/"
GAPs <- c("SpRNA1")
GAP.stock <- 10  # nM
GAP.working <- 1 # nM
Rans.stock.conc <- c(58)
sensor <- c('sensor')
sensor.stock.conc <- 125
reaction.volume <- 100 # ul
Ran.premix.volume <- 110
final.Ran.conc <- list()
#run_conc <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 12)
#final.Ran.conc[[1]] = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

#final.Ran.conc[[1]] = c(.25, .75, 1, 1.5, 2, 3)
#final.Ran.conc[[1]] = c(1, 2, 4, 6, 8, 12)
#final.Ran.conc[[1]] = c(1.5, 3, 3, 5, 7, 9, 11, 15)

#final.Ran.conc[[1]] = c(.25, .5, .75, 1, 1.25, 1.5, 2, 0, 4, 6, 10, 0)
#final.Ran.conc[[1]] = c(.25, 0.25, 0.5, 0.75, 0.75, 1, 1, 1.25, 1.5, 1.75, 2, 0, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 0)
final.Ran.conc[[1]] = c(16, 16, 16, 8, 8)





GAP.V <- round(reaction.volume/(GAP.stock/GAP.working), 2)
Ran.premix.to.add <- reaction.volume - GAP.V

reaction.mix.table <- data.frame()
index.table <- data.frame()
well_count <- 0
for (i in 1:length(Rans)) {
  prot <- Rans[i]
  for (j in seq_along(final.Ran.conc[[i]])) {
    well_count <- well_count + 1
    fin.Ran.conc <- final.Ran.conc[[i]][j]
    if (fin.Ran.conc < 3) {
      sensor.working <- 10 # uM
    } else if (fin.Ran.conc < 18) {
      sensor.working <- 20
    } else {
      sensor.working <- 30
    }
    Ran.premix.conc <- round(fin.Ran.conc * (1 + ( (reaction.volume - Ran.premix.to.add) / reaction.volume )), 2)
    stock.conc <- Rans.stock.conc[i]
    Ran.volume <- round(Ran.premix.volume/(stock.conc/Ran.premix.conc), 2)
    
    sensor.premix.conc <- round(sensor.working * (1 + ( (reaction.volume - Ran.premix.to.add) / reaction.volume )), 2)
    #sensor.V = reaction.volume/(sensor.stock.conc/sensor.working)
    sensor.V = round( Ran.premix.volume / (sensor.stock.conc / sensor.premix.conc), 2)
    
    
    
    Ran.premix.buffer.volume <- round(Ran.premix.volume - Ran.volume - sensor.V) 
    reaction.mix.table <- rbind(reaction.mix.table, data.frame(plate, "well" = paste0(row, well_count), prot, fin.Ran.conc, Ran.volume, 
                                                              sensor.V, Ran.premix.buffer.volume, 
                                                              GAP.stock, GAP.working, 
                                                              Ran.premix.to.add, GAP.V, date))
    index.table <- rbind(index.table, data.frame(plate, prot, "well" = paste0(row, well_count), fin.Ran.conc, GAP.working, date, "cutoff_time" = 3000))
  }
}
### (j+((i-1)*6))
reaction.mix.table
write.table(reaction.mix.table, file = paste0(outdir, "reaction_mix_temp.txt"), quote = F, row.names = F, sep = "\t")
(Ran.volume.sum <- reaction.mix.table %>%
    group_by(prot) %>% 
    summarise( "total_Ran_V" = sum(Ran.volume)) )
(sensor.volume.sum <- reaction.mix.table %>%
    group_by(prot) %>% 
    summarise( "total_diluted_125mM_sensor_V" = 1.05 * sum(sensor.V)) )   ### prepare 5 % extra diluted sensor
index.table
write.table(index.table, file = paste0(outdir, "index_table_temp.txt"), quote = F, row.names = F, sep = "\t")


