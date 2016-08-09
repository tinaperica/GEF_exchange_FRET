library(ggplot2)
library(minpack.lm)
setwd("~/Documents/GSP1_data/GEF_exchange_exp/")
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")
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
#proteins <- proteins[proteins != "PE7_Q147E"]
#proteins <- c("PE1_WT", "PE2_T34A")

get.initial.rate.subset <- function (total, conc) {
  min.total.fluor <- quantile(total$fluorescence, 0.001)
  max.total.fluor <- quantile(total$fluorescence, 0.999)
  total <- cbind(total, "percent.norm.fluorescence" = (max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor))
  total <- cbind(total, "substrate_conc" = conc - conc*((max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor)))
  cutoff.rxn.start.time <- numeric()
  if (conc < 1) {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.4])
  } else {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.5])
  }
  cutoff.rxn.start.percent <- subset(total, time < cutoff.rxn.start.time)
  return(cutoff.rxn.start.percent)
}
fit_rate_and_plot <- function (initial, prot, conc, GEF_conc) {
  time <- initial$time
  getPred <- function(pars, xx) pars$f0 * exp(- pars$k * time) + pars$c
  plot(initial$time, initial$substrate_conc, xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, k = 0.01, c = conc/2)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial$substrate_conc, xx = time, control = nls.lm.control(nprint = 1))
  lines(initial$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=2)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, substrate_conc = predicted.line, conc = conc)
  k <- coef(nls.out)[2]
  f0 <- coef(nls.out)[1]
  v0 = (f0 * k * exp(k * 0)) / (GEF_conc*0.001)
  return(data.frame(prot, conc, f0, k, v0, GEF_conc))
}

# fit_MM <- function (initial_vs_conc, prot) {
#   test.conc <- seq(min(initial_vs_conc$conc), max(initial_vs_conc$conc), 0.25)
#   conc <- initial_vs_conc$conc
#   getFit <- function(pars, xx) conc * pars$Vmax / (pars$Km + conc)
#   residFun <- function(p, observed, xx) observed - getFit(p,xx)
#   Km.estimate <- 1
#   Vm.estimate <- max(initial_vs_conc$v0)
#   parStart <- list(Vmax = Vm.estimate, Km = Km.estimate)
#   nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial_vs_conc$v0, xx = conc, control = nls.lm.control(nprint = 1))
#   getPred <- function(pars, xx) test.conc * pars$Vmax / (pars$Km + test.conc)
#   predicted.line <- getPred(as.list(coef(nls.out)), test.conc)
#   predicted.line.data.frame <- data.frame(conc = test.conc, v0 = predicted.line, prot = prot)
#   Vmax <- coef(nls.out)[1]
#   Km <- coef(nls.out)[2]
#   kcat = Vmax
#   parameters_and_predicted_line <- list("parameters" = data.frame(prot, Vmax, Km, kcat), "predicted_line_df" = predicted.line.data.frame)
#   return(parameters_and_predicted_line)
# }


fitting.parameters <- data.frame()
for (prot in proteins) {
  filename = paste(prot, outsufix, collapse = "_")
  pdf(file = filename)
  substr.concentrations <- unique(tecan$Gsp1_concentration[tecan$protein == prot])
  for (conc in substr.concentrations) {
    GEF_concentrations <- unique(tecan$GEF_concentration[tecan$protein == prot & tecan$Gsp1_concentration == conc])
    for (GEF_conc in GEF_concentrations) {
      data.subset <- subset(tecan, protein == prot & Gsp1_concentration == conc & GEF_concentration == GEF_conc, select=c(time = time, fluorescence = fluorescence))
      relevant.subset <- get.initial.rate.subset(data.subset, conc)
      fitting.parameters <- rbind(fitting.parameters, fit_rate_and_plot(relevant.subset, prot, conc, GEF_conc))
    }
  }
  dev.off()
}
fitting.parameters <- cbind(fitting.parameters,paste(fitting.parameters[[6]], "nM", sep = " "))
names(fitting.parameters)[length(names(fitting.parameters))] <- "GEF_conc_condition"

MM.parameters <- list()
KmVm.param.df <- data.frame()
KmVm.curve.to.plot <- data.frame()
Kmkcat.df <- data.frame()
for (i in 1:length(proteins)) {
  mutant = proteins[i]
  data.subset <- subset(fitting.parameters, prot == mutant, select=c(conc = conc, v0 = v0, GEF_conc = GEF_conc))
  MM.parameters[i] <- list (mutant = fit_MM(data.subset, mutant))
  KmVm.param.df <- rbind(KmVm.param.df, MM.parameters[[i]]$parameters)
  KmVm.curve.to.plot <- rbind(KmVm.curve.to.plot, MM.parameters[[i]]$predicted_line_df)
}
plot <- ggplot(data = fitting.parameters, aes(x = conc, y = v0, color = prot)) + geom_point()
plot <- plot + geom_line(data = KmVm.curve.to.plot)
plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~Delta~'F'~s^-1))
filename <- paste("MichaelisMenten_", Sys.Date(), ".pdf", sep = "")
pdf(file = filename, width = 10)
print(plot)
dev.off()

KmVm.param.df$kcat_over_Km <- KmVm.param.df$kcat / KmVm.param.df$Km
pairwise_filename <- paste("pairwise_MM_plots_", Sys.Date(), ".pdf", sep = "")
pdf(file = pairwise_filename, width = 10)
for (i in 1:length(proteins)) {
  protein = proteins[i]
  wt <- "PE1_WT"
  if (protein != wt) {
    temp.data <- subset(fitting.parameters, (prot == protein | prot == wt))
    temp.curve <- subset(KmVm.curve.to.plot, (prot == protein | prot == wt))
    plot <- ggplot(data = temp.data, aes(x = conc, y = v0, group = prot, color = prot))
    plot <- plot + geom_line(data = temp.curve)
    plot <- plot + geom_point(data = temp.data, aes(x = conc, y = v0, group = interaction(prot, GEF_conc_condition), color = prot, shape = GEF_conc_condition))
    plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~Delta~'F'~s^-1))
    print(plot)
  }
}
dev.off()

