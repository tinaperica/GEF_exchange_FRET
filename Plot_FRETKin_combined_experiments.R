library(ggplot2)
### In this code I add the filter file: tells which data to keep and which to discard (plot both into separate files to check)
## Data that is kept is the one where the reaction has been completed (reached plateau), and the ones where the reaction is not completed are discarded
max.na.rm <- function (x) {
  max(x, na.rm = T)
}
min.na.rm <- function (x) {
  min(x, na.rm = T)
}
setwd("~/Documents/GSP1_experimental_data/GEF_exchange_exp/")
name_prefix <- Sys.Date()  # date prefix for all the output file names
rawoutfilename = paste (name_prefix, "_FRET_kinetics_combined_raw.pdf", collapse = "")
#normoutfilename = paste(name_prefix, "_FRET_kinetics_norm.pdf", collapse = "")
scaledoutfilename = paste(name_prefix, "_FRET_kinetics_combined_scaled.pdf", collapse = "")
tecan<-read.delim("data/FRET_kinetics.txt", stringsAsFactors = F, head = F)
### the file with exp to discard was manually made based on raw plots
data.to.discard <- read.delim("data/FRET_kinetics_data_to_discard.txt", head = T) ### defined as date of exp, prot, prot conc, and GEF conc
#####
names(tecan) <- c("exp_date", "time", "row", "col", "protein", "condition", "Gsp1_concentration", "GEF_concentration", "fluorescence")
head(tecan)
tecan <- cbind(tecan,paste(tecan[[8]], "nM", sep = " "))
names(tecan)[length(names(tecan))] <- "GEF_conc_condition"
tecan <- tecan[complete.cases(tecan),]
tecan <- subset(tecan, GEF_concentration > 0)
tecan<-cbind(tecan, paste(tecan$exp_date, tecan$protein, tecan$Gsp1_concentration, tecan$GEF_concentration, sep = "_"))
names(tecan)[length(names(tecan))] <- "unique"
data.to.discard <- cbind(data.to.discard, paste(data.to.discard$exp_date, data.to.discard$protein, data.to.discard$Gsp1_concentration, data.to.discard$GEF_concentration, sep = "_"))
names(data.to.discard)[length(names(data.to.discard))] <- "unique"
unique_identifiers_to_discard <- as.character(data.to.discard$unique)
tecan.discarded <- tecan[tecan$unique %in% unique_identifiers_to_discard,]
tecan <- tecan[! tecan$unique %in% unique_identifiers_to_discard,]
proteins<-unique(tecan$protein)

tecan.start.max <- with(tecan, aggregate(fluorescence, by = list(protein = protein, condition = condition), max.na.rm))
tecan.start.min <- with(tecan, aggregate(fluorescence, by = list(protein = protein, condition = condition), min.na.rm))
tecan.normalized.max<-merge(tecan, tecan.start.max, by = c("protein", "condition"))
tecan.normalized <- merge(tecan.normalized.max, tecan.start.min, by = c("protein", "condition"))
names(tecan.normalized)[c(length(names(tecan.normalized)) - 1, length(names(tecan.normalized)))] <- c("max", "min") 
head(tecan.normalized)
tecan.normalized$normalized <- 1 - (tecan.normalized$max - tecan.normalized$fluorescence)/(tecan.normalized$max - tecan.normalized$min)
tecan.normalized$scaled <- (tecan.normalized$fluorescence - tecan.normalized$max)
tecan.normalized <- tecan.normalized[order(tecan.normalized$protein, tecan.normalized$Gsp1_concentration, tecan.normalized$time),]

plot_raw <- function (df, pref) {
  plots<-list()
  for (i in 1:length(proteins)) {
    prot <- proteins[i]
    temp<-subset(df, protein == prot)
    temp <- temp[order(temp$time),]
    plots[[i]] <- ggplot(data=temp, aes(x = time, y = fluorescence, group = interaction(condition, GEF_conc_condition), colour = condition, shape = GEF_conc_condition)) + geom_point(size = 0.7) #+ geom_line()
    plots[[i]] <- plots[[i]] + ggtitle(prot)
  }
  outfile = paste(pref, rawoutfilename, sep = "_") 
  pdf(file = outfile)
  bquiet<-lapply(plots, print)
  dev.off()
}
plot_raw (tecan, "data_to_keep")
plot_raw(tecan.discarded, "discarded_data")


# plots<-list()
# for (i in 1:length(proteins)) {
#   prot <- proteins[i]
#   temp<-subset(tecan.normalized, protein == prot)
#   temp <- temp[order(temp$time),]
#   plots[[i]]<-ggplot(data=temp, aes(x = time, y = normalized, group = interaction(condition, GEF_conc_condition), colour = condition, shape = GEF_conc_condition)) + geom_point()
#   #plots[[i]]<-plots[[i]] + stat_smooth(span = 0.5)
#   plots[[i]]<-plots[[i]] + ggtitle(prot)
# }
# pdf(file = normoutfilename)
# bquiet<-lapply(plots, print)
# dev.off()
plots<-list()
for (i in 1:length(proteins)) {
  prot <- proteins[i]
  temp<-subset(tecan.normalized, protein == prot)
  temp <- temp[order(temp$time),]
  plots[[i]]<-ggplot(data=temp, aes(x = time, y = scaled, group = interaction(condition, GEF_conc_condition), colour = condition, shape = GEF_conc_condition)) + geom_point()
  #plots[[i]]<-plots[[i]] + stat_smooth(span = 0.5)
  plots[[i]]<-plots[[i]] + ggtitle(prot)
}
pdf(file = scaledoutfilename)
bquiet<-lapply(plots, print)
dev.off()
