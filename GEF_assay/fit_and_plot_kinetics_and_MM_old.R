library(ggplot2)
library(minpack.lm)
setwd("~/Documents/GSP1_experimental_data/GEF_exchange_exp/")
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")
tecan<-read.delim("data/FRET_kinetics.txt", stringsAsFactors = F, head = F)
### the file with exp to discard was manually made based on inspection of all the raw data plots
data.to.discard <- read.delim("data/FRET_kinetics_data_to_discard.txt", head = T) ### defined as date of exp, prot, prot conc, and GEF conc
#####
names(tecan) <- c("exp_date", "time", "row", "col", "protein", "condition", "Gsp1_concentration", "GEF_concentration", "fluorescence")
head(tecan)
tecan <- tecan[tecan$Gsp1_concentration < 10,]
tecan <- cbind(tecan,paste(tecan[[8]], "nM", sep = " "))
names(tecan)[length(names(tecan))] <- "GEF_conc_condition"
tecan <- tecan[complete.cases(tecan),]
tecan <- subset(tecan, GEF_concentration > 0)
tecan<-cbind(tecan, paste(tecan$exp_date, tecan$protein, tecan$Gsp1_concentration, tecan$GEF_concentration, sep = "_"))
names(tecan)[length(names(tecan))] <- "unique"  #### unique is a combination of experiment date and conditions (Ran and GEF concentrations)
data.to.discard <- cbind(data.to.discard, paste(data.to.discard$exp_date, data.to.discard$protein, data.to.discard$Gsp1_concentration, data.to.discard$GEF_concentration, sep = "_"))
names(data.to.discard)[length(names(data.to.discard))] <- "unique"
unique_identifiers_to_discard <- as.character(data.to.discard$unique)
tecan.discarded <- tecan[tecan$unique %in% unique_identifiers_to_discard,]
tecan <- tecan[! tecan$unique %in% unique_identifiers_to_discard,]
proteins<-unique(tecan$protein)

proteins <- c("PE1_WT", "PE3_T34E", "PE2_T34A", "PE9_D79S", "PE10_G80A", "PE4_R108L", "PE18_R108Y", "PE7_Q147E", "PE14_Y157A", "PE29_R112S" )


get.initial.rate.subset <- function (total, conc, GEF_conc) {
  min.total.fluor <- quantile(total$fluorescence, 0.001)
  max.total.fluor <- quantile(total$fluorescence, 0.999)
  total <- cbind(total, "percent.norm.fluorescence" = (max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor))
  total <- cbind(total, "substrate_conc" = conc - conc*((max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor)))
  cutoff.rxn.start.time <- numeric()
  if (conc < 1 | GEF_conc > 30) {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.5]) #0.5
  } else {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.7])  # 0.7
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
  plot(initial$time, initial$substrate_conc, xlim = c(0, 8000), xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, exp, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, k = 0.01, c = conc/2)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial$substrate_conc, xx = time, control = nls.lm.control(nprint = 1))
  lines(initial$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, substrate_conc = predicted.line, conc = conc)
  k <- coef(nls.out)[2]
  f0 <- coef(nls.out)[1]
  v0 = (f0 * k * exp(k * 0)) / (GEF_conc*0.001)
  return(data.frame(prot, conc, f0, k, v0, GEF_conc, exp))
  par(opar)
}

fit_MM <- function (initial_vs_conc, prot) {
  test.conc <- seq(min(initial_vs_conc$conc), max(initial_vs_conc$conc), 0.25)
  conc <- initial_vs_conc$conc
  getFit <- function(pars, xx) conc * pars$Vmax / (pars$Km + conc)
  residFun <- function(p, observed, xx) observed - getFit(p,xx)
  Km.estimate <- 1
  Vm.estimate <- max(initial_vs_conc$v0)
  parStart <- list(Vmax = Vm.estimate, Km = Km.estimate)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial_vs_conc$v0, xx = conc, control = nls.lm.control(nprint = 1))
  getPred <- function(pars, xx) test.conc * pars$Vmax / (pars$Km + test.conc)
  predicted.line <- getPred(as.list(coef(nls.out)), test.conc)
  predicted.line.data.frame <- data.frame(conc = test.conc, v0 = predicted.line, prot = prot)
  Vmax <- round(coef(nls.out)[1], 3)
  Km <- round(coef(nls.out)[2], 3)
  kcat = Vmax
  
  standard_errors_of_parameters <- round(summary(nls.out)$coefficients[, 2], 3)
  kcat_error <- round(as.numeric(standard_errors_of_parameters[[1]]), 3)
  Km_error <- round(as.numeric(standard_errors_of_parameters[[2]]), 3)
  summary <- summary(nls.out)
  sd_of_fit <- round(as.numeric(summary[["sigma"]]), 3)
  parameters_and_predicted_line <- list("parameters" = data.frame(prot, Vmax, Km, kcat, kcat_error, Km_error, sd_of_fit), 
                                        "predicted_line_df" = predicted.line.data.frame)
  return(parameters_and_predicted_line)
}

fitting.parameters <- data.frame()
for (i in seq_along(proteins)) {
  prot <- proteins[i]
  filename = paste0(Sys.Date(), "_", prot, outsufix)
  pdf(file = filename, height = 10)
  substr.concentrations <- unique(tecan$Gsp1_concentration[tecan$protein == prot])
  for (j in seq_along(substr.concentrations))  {
    conc = substr.concentrations[j]
    GEF_concentrations <- unique(tecan$GEF_concentration[tecan$protein == prot & tecan$Gsp1_concentration == conc])
    for (k in seq_along(GEF_concentrations)) {
      GEF_conc = GEF_concentrations[k]
      experiments <- unique(tecan$exp_date[tecan$protein == prot & tecan$Gsp1_concentration == conc & tecan$GEF_concentration == GEF_conc])
      for (l in seq_along(experiments)) {
        exp <- experiments[l]
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
#plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~Delta~'F'~s^-1))
plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~s^-1))
filename <- paste0("MichaelisMenten_", Sys.Date(), ".pdf")
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
    #plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~Delta~'F'~s^-1))
    plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~s^-1))
    print(plot)
  }
}
dev.off()





# fit_kobs <- function (kobs_vs_conc, prot) {
#   test.conc <- seq(min(kobs_vs_conc), max(kobs_vs_conc), 0.25)
#   conc <- kobs_vs_conc$conc
#   GEF_conc <- kobs_vs_conc$GEF_conc
#   getFit <- function (pars, xx) pars$knuc + ( (pars$K3*pars$kdiss*GEF_conc) / (1 + pars$K3*conc) )
#   residFun <- function (p, observed, xx) observed - getFit(p,xx)
#   parStart <- list (knuc = 0.01, K3 = 1, kdiss = 5)
#   nls.out <- nls.lm(par = parStart, fn = residFun, observed = kobs_vs_conc$k, xx = conc, control = nls.lm.control(nprint = 1))
#   getPred <- function(pars, xx) pars$knuc + ( (pars$K3*pars$kdiss*GEF_conc) / (1 + pars$K3*test.conc) )
#   predicted.line <- getPred (as.list(coef(nls.out)), test.conc)
#   predicted.line.data.frame <- data.frame(conc = test.conc, kobs = predicted.line, prot = prot)
#   knuc <- coef(nls.out)[1]
#   K3 <- coef(nls.out)[2]
#   kdiss <- coef(nls.out)[3]
#   parameters_and_predicted_line <- list("parameters" = data.frame(prot, knuc, K3, kdiss), "predicted_line_df" = predicted.line.data.frame)
#   return(parameters_and_predicted_line)
# }
# kobs.parameters <- list()
# kobs.param.df <- data.frame()
# kobs.curve.to.plot <- data.frame()
# kobs.df <- data.frame()
# for (i in 1:length(proteins)) {
#   mutant = proteins[i]
#   data.subset <- subset(fitting.parameters, prot == mutant, select=c(conc = conc, k = k, GEF_conc = GEF_conc, exp = exp))
#   kobs.parameters[i] <- list (mutant = fit_kobs(data.subset, mutant))
#   kobs.param.df <- rbind(kobs.param.df, kobs.parameters[[i]]$parameters)
#   kobs.curve.to.plot <- rbind(kobs.curve.to.plot, kobs.parameters[[i]]$predicted_line_df)
# }
# plot <- ggplot(data = fitting.parameters, aes(x = conc, y = k, color = prot)) + geom_point()
# plot <- plot + geom_line(data = kobs.curve.to.plot)
# plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('kobs / ~s^-1'))
# filename <- paste("kobs_", Sys.Date(), ".pdf", sep = "")
# pdf(file = filename, width = 10)
# print(plot)
# dev.off()


Km.barplot <- barplot(height = KmVm.param.df$Km, names.arg = KmVm.param.df$prot, 
                      main = "Km of GEF exchange", ylab = bquote("Km / " *mu~M),
                      xaxt = "n", cex.names = 0.75, ylim = c(0,13))
text(x = Km.barplot, y = par("usr")[3] - 1, srt = 45, adj = 1, labels = KmVm.param.df$prot, xpd = TRUE)
segments(Km.barplot, KmVm.param.df$Km - KmVm.param.df$Km_error, Km.barplot,
          KmVm.param.df$Km + KmVm.param.df$Km_error, lwd = 1.5)
arrows(Km.barplot, KmVm.param.df$Km - KmVm.param.df$Km_error, Km.barplot,
       KmVm.param.df$Km + KmVm.param.df$Km_error, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)



kcat.barplot <- barplot(height = KmVm.param.df$kcat, names.arg = KmVm.param.df$prot, 
                      main = "kcat of GEF exchange", ylab = bquote("kcat / s-1"),
                      xaxt = "n", cex.names = 0.6, ylim = c(0,4))
text(x = kcat.barplot, y = par("usr")[3] - 0.3, srt = 45, adj = 1, labels = KmVm.param.df$prot, xpd = TRUE)
segments(kcat.barplot, KmVm.param.df$kcat - KmVm.param.df$kcat_error, kcat.barplot,
         KmVm.param.df$kcat + KmVm.param.df$kcat_error, lwd = 1.5)
arrows(kcat.barplot, KmVm.param.df$kcat - KmVm.param.df$kcat_error, kcat.barplot,
       KmVm.param.df$kcat + KmVm.param.df$kcat_error, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)





