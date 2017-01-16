library(ggplot2)
max.na.rm <- function (x) {
  max(x, na.rm = T)
}
min.na.rm <- function (x) {
  min(x, na.rm = T)
}
setwd("~/Documents/GSP1_experimental_data/GEF_exchange_exp/")
name_prefix <- "20160830"
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")
tecan<-read.delim("data/20160830_FRET_kinetics_parsed.txt", stringsAsFactors = F, head = F)
names(tecan) <- c("exp_date", "time", "row", "col", "protein", "condition", "Gsp1_concentration", "GEF_concentration", "fluorescence")
tecan <- cbind(tecan,paste(tecan[[8]], "nM", sep = " "))
names(tecan)[length(names(tecan))] <- "GEF_conc_condition"
tecan <- subset(tecan, time < 4000)
tecan <- tecan[complete.cases(tecan),]
rawoutfilename = paste (name_prefix, "_FRET_kinetics_raw.pdf", collapse = "")
scaledoutfilename = paste(name_prefix, "_FRET_kinetics_scaled.pdf", collapse = "")
proteins<-unique(tecan$protein)
#### define fitting functions
get.initial.rate.subset <- function (total, conc, GEF_conc) {
  min.total.fluor <- quantile(total$fluorescence, 0.001)
  max.total.fluor <- quantile(total$fluorescence, 0.999)
  total <- cbind(total, "percent.norm.fluorescence" = (max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor))
  total <- cbind(total, "substrate_conc" = conc - conc*((max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor)))
  #### these are empirical thresholds - check the actual fitting files and modify accordingly if necessary
  cutoff.rxn.start.time <- numeric()
  if ( GEF_conc < 15) {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.7]) 
  } else if (GEF_conc >= 15 & GEF_conc < 30) {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.6])  
  } else {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.6])
  }
  cutoff.rxn.start.percent <- subset(total, time < cutoff.rxn.start.time)
  return(cutoff.rxn.start.percent)
}
fit_rate_and_plot <- function (all, initial, prot, conc, GEF_conc, exp) {
  time <- initial$time
  getPred <- function(pars, xx) pars$f0 * exp(- pars$k * time) + pars$c
  opar <- par()
  op<-par(mfrow=c(2,1))
  plot(all$time, all$fluorescence, xlab = "time / s", ylab = "fluorescence", main = paste(prot, conc, GEF_conc, exp, collapse = ""))
  plot(initial$time, initial$substrate_conc, xlim = c(0, 2000), xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, exp, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, k = 0.01, c = conc/2)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial$substrate_conc, xx = time, control = nls.lm.control(nprint = 1))
  lines(initial$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, substrate_conc = predicted.line, conc = conc)
  k <- coef(nls.out)[2]
  f0 <- coef(nls.out)[1]
  v0 = (f0 * k * exp(k * 0))  #### for Michaelis Menten this v0 is immediately divided by the GEF concentration
  return(data.frame(prot, conc, f0, k, v0, GEF_conc, exp))
  par(opar)
}

tecan.start.max <- with(tecan, aggregate(fluorescence, by = list(protein = protein, GEF_concentration = GEF_concentration), max.na.rm))
tecan.start.min <- with(tecan, aggregate(fluorescence, by = list(protein = protein, GEF_concentration = GEF_concentration), min.na.rm))
temp.max<-merge(tecan, tecan.start.max, by = c("protein", "GEF_concentration"))
tecan.scaled <- merge(tecan.start.min, temp.max, by = c("protein", "GEF_concentration"))
names(tecan.scaled)[c(3, 12)] <- c("min", "max") 
tecan.scaled$scaled <- tecan.scaled$fluorescence - tecan.normalized$max
#### this minmaxpercent is used to determine which conditions are just photobleaching
## because the exchange has already finished before the measurement starts (this is with most high conc of GEF with WT Gsp1)
tecan.scaled$minmaxpercent <- (tecan.scaled$max - tecan.scaled$min) / tecan.scaled$max
plot(density(tecan.scaled$minmaxpercent[tecan.scaled$GEF_concentration == 0]))
plot(density(tecan.scaled$minmaxpercent[tecan.scaled$GEF_concentration > 0]))
tecan.zero <- subset(tecan.scaled, GEF_concentration == 0)
tecan.clean <- data.frame()
for (prot in proteins) {
  temp <- subset(tecan.scaled, protein == prot)
  refminmaxpercent <- tecan.zero$minmaxpercent[tecan.zero$protein == prot]
  refminmaxpercent <- 2*refminmaxpercent
  relevant.tecan.temp <- subset(temp, minmaxpercent > refminmaxpercent)
  tecan.clean <- rbind(tecan.clean, relevant.tecan.temp)
}
tecan.clean <- rbind(tecan.clean, tecan.zero)
tecan.clean <- tecan.clean[order(tecan.clean$protein, tecan.clean$GEF_concentration, tecan.clean$time),]
tecan <- tecan.clean
plots<-list()
for (i in 1:length(proteins)) {
  prot <- proteins[i]
  temp<-subset(tecan, protein == prot)
  temp <- temp[order(temp$time),]
  plots[[i]] <- ggplot(data=temp, aes(x = time, y = fluorescence, group = interaction(condition, GEF_conc_condition), colour = GEF_conc_condition, shape = condition)) + geom_point(size = 0.7) #+ geom_line()
  plots[[i]] <- plots[[i]] + ggtitle(prot)
}
pdf(file = rawoutfilename)
bquiet<-lapply(plots, print)
dev.off()

plots<-list()
for (i in 1:length(proteins)) {
  prot <- proteins[i]
  temp<-subset(tecan.normalized, protein == prot)
  temp <- temp[order(temp$time),]
  plots[[i]]<-ggplot(data=temp, aes(x = time, y = scaled, group = interaction(condition, GEF_conc_condition), colour = GEF_conc_condition, shape = condition)) + geom_point()
  #plots[[i]]<-plots[[i]] + stat_smooth(span = 0.5)
  plots[[i]]<-plots[[i]] + ggtitle(prot)
}
pdf(file = scaledoutfilename)
bquiet<-lapply(plots, print)
dev.off()

fitting.parameters <- data.frame()
for (prot in proteins) {
  filename = paste(prot, outsufix, collapse = "_")
  pdf(file = filename, height = 10)
  substr.concentrations <- unique(tecan$Gsp1_concentration[tecan$protein == prot])
  for (conc in substr.concentrations) {
    GEF_concentrations <- unique(tecan$GEF_concentration[tecan$protein == prot & tecan$Gsp1_concentration == conc])
    for (GEF_conc in GEF_concentrations) {
      experiments <- unique(tecan$exp_date[tecan$protein == prot & tecan$Gsp1_concentration == conc & tecan$GEF_concentration == GEF_conc])
      for (exp in experiments) {
        data.subset <- subset(tecan, protein == prot & Gsp1_concentration == conc & GEF_concentration == GEF_conc & tecan$exp_date == exp, select=c(time = time, fluorescence = fluorescence))
        relevant.subset <- get.initial.rate.subset(data.subset, conc, GEF_conc)
        fitting.parameters <- rbind(fitting.parameters, fit_rate_and_plot(data.subset, relevant.subset, prot, conc, GEF_conc, exp))
      }
    }
  }
  dev.off()
}
fitting.parameters <- cbind(fitting.parameters,paste(fitting.parameters[[6]], "nM", sep = " "))
names(fitting.parameters)[length(names(fitting.parameters))] <- "GEF_conc_condition"
fitting.parameters <- fitting.parameters[order(fitting.parameters$prot, fitting.parameters$GEF_conc),]
colors <- rainbow(length(proteins))
plot(fitting.parameters$GEF_conc, fitting.parameters$v0, xlab = "SRM1 concentration / nM", ylab = "v0 / [Gsp1:GDP] uM s-1", type = "n")
for (i in 1:length(proteins)) {
  prot = proteins[i]
  temp.subset <- fitting.parameters[fitting.parameters$prot == prot,]
  lines(temp.subset$GEF_conc, temp.subset$v0, type = "o", col = colors[i])
}
legend("top", fill = colors, legend = proteins, bty = "n")



