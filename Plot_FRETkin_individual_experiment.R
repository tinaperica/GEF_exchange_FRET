library(ggplot2)
max.na.rm <- function (x) {
  max(x, na.rm = T)
}
min.na.rm <- function (x) {
  min(x, na.rm = T)
}
setwd("~/Documents/GSP1_experimental_data/GEF_exchange_exp/")
name_prefix <- "20160906"
tecan<-read.delim("data/20160906_FRET_kinetics_parsed.txt", stringsAsFactors = F, head = F)
names(tecan) <- c("exp_date", "time", "row", "col", "protein", "condition", "Gsp1_concentration", "GEF_concentration", "fluorescence")
head(tecan)
tecan <- cbind(tecan,paste(tecan[[8]], "nM", sep = " "))
names(tecan)[length(names(tecan))] <- "GEF_conc_condition"
tecan <- tecan[complete.cases(tecan),]
tecan <- subset(tecan, GEF_concentration > 0)
rawoutfilename = paste (name_prefix, "_FRET_kinetics_raw.pdf", collapse = "")
#normoutfilename = paste(name_prefix, "_FRET_kinetics_norm.pdf", collapse = "")
scaledoutfilename = paste(name_prefix, "_FRET_kinetics_scaled.pdf", collapse = "")
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


plots<-list()
for (i in 1:length(proteins)) {
  prot <- proteins[i]
  temp<-subset(tecan, protein == prot)
  temp <- temp[order(temp$time),]
  plots[[i]] <- ggplot(data=temp, aes(x = time, y = fluorescence, group = interaction(condition, GEF_conc_condition), colour = condition, shape = GEF_conc_condition)) + geom_point(size = 0.7) #+ geom_line()
  plots[[i]] <- plots[[i]] + ggtitle(prot)
  # for (j in 1:length(unique(temp$Gsp1_concentration))) {
  #   conc <- unique(temp$Gsp1_concentration)[j]
  #   temp.conc <- subset(temp, Gsp1_concentration == conc)
  #   temp.conc <- temp.conc[order(temp.conc$time),]
  # }
}
pdf(file = rawoutfilename)
bquiet<-lapply(plots, print)
dev.off()
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

