library(tidyverse)
SubstrateDilute <- function(stock, work, volume, Mw) {
  molar.stock <- 1000*(stock/Mw)
  dilutionX <- molar.stock/work
  protein.V <- round(volume/dilutionX,1)
  buffer.V <- round(volume - protein.V,1)
  return(data.frame(protein.V, buffer.V, "dilution" = round(dilutionX, 2)))
}

### Set up  experiment parameters
plate <- 4
row <- "A"
date <- "20180213"
outdir <- "GEF_assay/biotek_test_exp/"
Rans <- c("PE9_D79S", "PE10_G80A")
GEFs <- c("SRM1")
GEF.stock <- 500  # nM
GEF.working <- 5 # nM
Rans.stock.conc <- c(40, 40)
nucleotides<-c('mant-GTP')
nucleotide.stock.conc <- 2000
#mant.stock <- 2000 # uM
nucleotide.working <- 200 # uM
reaction.volume <- 100 # ul
Ran.premix.volume <- 110
final.Ran.conc <- list()
final.Ran.conc[[1]] = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 9, 10, 11)
final.Ran.conc[[2]] = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 9, 10, 11)
#final.Ran.conc[[1]] = c(1.75, 2, 2.5, 3, 3.75, 4)
#final.Ran.conc[[2]] = c(1.75, 2, 2.5, 3, 3.75, 4)
#final.Ran.conc[[1]] = c(0.25, 0.75, 1, 1.5, 2.25, 3.25)
#final.Ran.conc[[2]] = c(0.25, 0.75, 1, 1.5, 2.25, 3.25)
#final.Ran.conc[[1]] = c(0.5, 1.25, 2.75, 3.5, 4.5, 5)
#final.Ran.conc[[2]] = c(0.5, 1.25, 2.75, 3.5, 4.5, 5)
#final.Ran.conc[[1]] = c(0.5, 1, 1.5, 2, 2.5, 3, 5, 6, 8, 10, 11, 13)

nucleotide.V = reaction.volume/(nucleotide.stock.conc/nucleotide.working)
GEF.V <- round(reaction.volume/(GEF.stock/GEF.working), 2)
GEF.nucleotide.premix.volume <- GEF.V + nucleotide.V
Ran.premix.to.add <- reaction.volume - GEF.nucleotide.premix.volume

reaction.mix.table <- data.frame()
index.table <- data.frame()
for (i in 1:length(Rans)) {
  prot <- Rans[i]
  for (j in seq_along(final.Ran.conc[[i]])) {
    fin.Ran.conc <- final.Ran.conc[[i]][j]
    Ran.premix.conc <- round(fin.Ran.conc * (1 + ( (reaction.volume - Ran.premix.to.add) / reaction.volume )), 2)
    stock.conc <- Rans.stock.conc[i]
    Ran.volume <- round(Ran.premix.volume/(stock.conc/Ran.premix.conc), 2)
    Ran.premix.buffer.volume <- round(Ran.premix.volume - Ran.volume) 
    reaction.mix.table <- rbind(reaction.mix.table, data.frame(plate, "well" = paste0(row,(j+((i-1)*6))), prot, fin.Ran.conc, Ran.volume, Ran.volume/2, Ran.premix.buffer.volume, 
                                                              GEF.stock, GEF.working, 
                                                              Ran.premix.to.add, GEF.nucleotide.premix.volume, date))
    index.table <- rbind(index.table, data.frame(plate, prot, "well" = paste0(row, (j+((i-1)*6))), fin.Ran.conc, GEF.working, date, "cutoff_time" = 3000))
  }
}

reaction.mix.table
write.table(reaction.mix.table, file = paste0(outdir, "reaction_mix_temp.txt"), quote = F, row.names = F, sep = "\t")
(Ran.volume.sum <- reaction.mix.table %>%
    group_by(prot) %>% 
    summarise( "total_Ran_V" = sum(Ran.volume)) )

index.table
write.table(index.table, file = paste0(outdir, "index_table_temp.txt"), quote = F, row.names = F, sep = "\t")


