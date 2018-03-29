#setwd("~/Documents/GSP1_data/GEF_exchange_exp/")
SubstrateDilute <- function(stock, work, volume, Mw) {
  molar.stock <- 1000*(stock/Mw)
  dilutionX <- molar.stock/work
  protein.V <- round(volume/dilutionX,1)
  buffer.V <- round(volume - protein.V,1)
  return(data.frame(protein.V, buffer.V, "dilution" = round(dilutionX, 2)))
}
frozen.stock.table <- read.delim("current_protein_stocks_in_mg_per_ml.txt", head = F)
names(frozen.stock.table) <- expression(number, mutant, concentration)
frozen.stock.table[["protein"]] <- paste0("PE", frozen.stock.table$number, "_", frozen.stock.table$mutant)
setwd("~/Documents/GSP1/GSP1_experimental_data/GEF_exchange_exp")
### Set up  experiment parameters
GEF.stock <- list()
GEF.working <- list()
GEFs <- c("SRM1")
GEF.stock[[1]] <- c(10)  # nM
GEF.working[[1]] <- c(5) # nM
#GEF.stock[[2]] <- c(500)  # nM
#GEF.working[[2]] <- c(250) # nM
#GEF.stock[[3]] <- c(0)
#GEF.working[[3]] <-c(0)
Rans <- c("PE1_WT")
Rans.stock.conc <- c(39.3)
nucleotides<-c('mant-dGDP')
nucleotide.stock <- c(2000)
#mant.stock <- 2000 # uM
nucleotide.working <- list()
nucleotide.working[[1]] <- c(200) # uM
# nucleotide.working[[2]] <- 200
# nucleotide.working[[3]] <- 200
# nucleotide.working[[4]] <- 200
# nucleotide.working[[5]] <- 0
reaction.volume <- 100 # ul
final.Ran.conc <- list()
final.Ran.conc[[1]] = c(0.5, 0.5, 2, 4, 7, 12, 12)

reaction.mix.table <- data.frame()
for (i in 1:length(Rans)) {
  prot <- Rans[i]
  for (j in 1:length(final.Ran.conc[[i]])) {
    Ran.conc <- final.Ran.conc[[i]][j]
    stock.conc <- Rans.stock.conc[i]
    Ran.volume <- round(reaction.volume/(stock.conc/Ran.conc), 2)
    for (g in 1:length(GEFs)) {
      for (k in 1:length(GEF.stock[[g]])) {
        for (n in 1:length(nucleotides)) {
          for (nw in 1:length(nucleotide.working[[n]])) {
            nucleotide.name <- nucleotides[n]
            nucleotide.stock.conc <- nucleotide.stock[n]
            if (nucleotide.stock.conc != 0) {
              nucleotide.V = reaction.volume/(nucleotide.stock.conc/nucleotide.working[[n]][nw])
            } else {
              nucleotide.V <- 0
              
            }
            if (GEF.stock[[g]][k] != 0) {
              GEF.V <- round(reaction.volume/(GEF.stock[[g]][k]/GEF.working[[g]][k]), 2)
            } else {
              GEF.V <- 0
            }
            buffer.volume = round(reaction.volume - (Ran.volume + nucleotide.V + GEF.V), 2)
            GEF.stock.name <- as.character(GEF.stock[[g]][k])
            reaction.mix.table <- rbind(reaction.mix.table, data.frame(prot, Ran.conc, Ran.volume, Ran.volume/2, buffer.volume, buffer.volume/2, nucleotide.name, nucleotide.V, GEF.stock.name, "GEF.conc" = round(GEF.working[[g]][k],2), GEF.V))
          }
        }
      }
    }
  }
}

reaction.mix.table
prot.stocks.table <- data.frame()
for (i in 1:length(Rans)) {
  stock <- sum(reaction.mix.table$Ran.volume[reaction.mix.table$prot == Rans[i]])
  prot.stocks.table <- rbind(prot.stocks.table, data.frame("prot" = Rans[i], "stock.V" = stock))
}
sum(reaction.mix.table$buffer.volume)
sum(reaction.mix.table$GEF.V)
sum(reaction.mix.table$nucleotide.V)
prot.stocks.table

#### Make Stocks
MakeStocks.df <- data.frame()
for (i in seq_along(Rans)) {
  prot <- Rans[i]
  frozen.stock.conc <- frozen.stock.table$conc[frozen.stock.table$prot == prot]
  Stock <- prot.stocks.table$stock.V[prot.stocks.table$prot == prot]
  StockToMake <- Stock + 0.5*Stock  ### make 50% more
  dilution.table <- SubstrateDilute(frozen.stock.conc, 35, StockToMake, 24.8)
  MakeStocks.df <- rbind(MakeStocks.df, data.frame("prot" = prot, dilution.table))
}
MakeStocks.df
### make a serial dilution table for the enzyme, given the uM concentration of the original enzyme stock,
### a range of enzyme nM concentrations (from the GEF.stock list) -> this range WILL be sorted in descending order
### and the total volume of enzyme (the function takes into account that part of the volume is used in the serial dilution)
Enzyme_Serial_Dilution <- function(conc_range_nM, original_stock_uM, EnzymeTotalVolume) {
  conc_range_nM <- sort(conc_range_nM, decreasing = T)
  dilution_names <- c("stock",LETTERS[1:length(conc_range_nM)])
  dilutions <- vector()
  dilution_table <- data.frame()
  dilutions <- append(dilutions, original_stock_uM*1000/conc_range_nM[1])
  for (i in 2:length(conc_range_nM)) {
    dilutions <- append(dilutions, round(conc_range_nM[i-1]/conc_range_nM[i],1))
  }
  passed_volume <- rep(0, length(conc_range_nM))
  for (i in length(conc_range_nM):2) {
    volume <- round(EnzymeTotalVolume + 0.2*EnzymeTotalVolume + passed_volume[i+1], 2)  ## passed_volume[i] is the volume passed on from the previous dilution
    passed_volume[i+1] <- round(volume/dilutions[i],2)  ### this is the volume passed to the next dilution
    dilution_table <- rbind(dilution_table, data.frame("name" = dilution_names[i], "fold dilution" = dilutions[i], "from sample" = dilution_names[i-1], "protein volume" = passed_volume[i-1], "buffer volume" = round(volume - passed_volume[i-1], 2), "final volume" = volume - passed_volume[i]))
  }
  dilution_table<-dilution_table[order(dilution_table$name, decreasing = T),]
}

### this function calculates the dilution of substrate protein with buffer given the stock mass concentration in mg/mol, the final molar concentration in uM, and Mw in kDa


