library(ggplot2)
library(minpack.lm)
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")
data<-read.delim("GEF_assay/data/20171018_TP_FRET_kinetics_parsed.txt", stringsAsFactors = F, head = F)
### the file with exp to discard was manually made based on inspection of all the raw data plots
#####
names(data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence)
head(data)

data$well <- do.call(paste0, data[c("row", "column")])

data.to.discard <- read.delim("GEF_assay/data/FRET_kinetics_data_to_discard.txt", head = T) ### defined as date of exp, prot, prot conc, and GEF conc
#####
data<-cbind(data, paste(data$date, data$sample, data$conc, data$GEF_conc, data$well, sep = " "))
names(data)[length(names(data))] <- "unique"  #### unique is a combination of experiment date and conditions (Ran and GAP conc and well)

df_args <- c(data.to.discard, sep=" ")
data.to.discard$unique <- do.call(paste, df_args)
#data <- cbind(data, paste(data[[8]], "nM", sep = " "))
#names(data)[length(names(data))] <- "GEF_conc_condition"
#tecan <- tecan[complete.cases(tecan),]
#tecan <- subset(tecan, GEF_conc > 0)
#tecan<-cbind(tecan, paste(tecan$exp_date, tecan$protein, tecan$conc, tecan$GEF_conc, sep = "_"))
#names(tecan)[length(names(tecan))] <- "unique"  #### unique is a combination of experiment date and conditions (Ran and GEF concentrations)
data.to.discard <- cbind(data.to.discard, paste(data.to.discard$exp_date, data.to.discard$protein, data.to.discard$conc, data.to.discard$GEF_conc, sep = "_"))
names(data.to.discard)[length(names(data.to.discard))] <- "unique"
unique_identifiers_to_discard <- as.character(data.to.discard$unique)
discarded.data <- data[data$unique %in% unique_identifiers_to_discard,]
data <- data[! data$unique %in% unique_identifiers_to_discard,]
proteins<-unique(data$sample)

get.initial.rate.subset <- function (total, conc, GEF_conc) {
  min.total.fluor <- quantile(total$fluorescence, 0.001)
  max.total.fluor <- quantile(total$fluorescence, 0.999)
  total <- cbind(total, "percent.norm.fluorescence" = (max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor))
  total <- cbind(total, "substrate_conc" = conc - conc*((max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor)))
  cutoff.rxn.start.time <- numeric()
  if (conc < 2) {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.4]) #0.5
  } else {
    cutoff.rxn.start.time <- max(total$time[total$percent.norm.fluorescence < 0.8])  # 0.7
  }
  cutoff.rxn.start.percent <- subset(total, time < cutoff.rxn.start.time)
  return(cutoff.rxn.start.percent)
}


get.substrate.conc <- function (total, conc, GEF_conc) {
  min.total.fluor <- quantile(total$fluorescence, 0.001)
  max.total.fluor <- quantile(total$fluorescence, 0.999)
  total <- cbind(total, "percent.norm.fluorescence" = (max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor))
  total <- cbind(total, "substrate_conc" = conc - conc*((max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor)))
  return(total)
}

fit_exp_decay_and_plot <- function (all, initial, prot, conc, GEF_conc, exp) {
  time <- initial$time
  getPred <- function(pars, xx) (pars$f0 - pars$c) * exp(- pars$k * time) + pars$c
  opar <- par()
  op<-par(mfrow=c(2,1))
  plot(all$time, all$fluorescence, xlab = "time / s", ylab = "fluorescence", main = paste(prot, conc, GEF_conc, exp, collapse = ""))
  plot(initial$time, initial$substrate_conc, xlim = c(0, 10000), xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, exp, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, k = 0.01, c = 0)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial$substrate_conc, xx = time, control = nls.lm.control(nprint = 1))
  lines(initial$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, substrate_conc = predicted.line, conc = conc)
  k <- coef(nls.out)[2]
  f0 <- coef(nls.out)[1]
  v0 = (f0 * k * exp(k * 0)) #/ (GEF_conc*0.001)
  return(data.frame(prot, conc, f0, k, v0, GEF_conc, exp))
  par(opar)
}

fit_two_phase_decay_and_plot <- function (all, initial, prot, conc, GEF_conc, exp) {
  time <- initial$time
  #SpanFast=(Y0-Plateau)*PercentFast*.01
  #SpanSlow=(Y0-Plateau)*(100-PercentFast)*.01
  #Y=Plateau + SpanFast*exp(-KFast*X) + SpanSlow*exp(-KSlow*X)
  
  #getPred <- function(pars, xx) pars$f0 * exp(- pars$k * time) + pars$fmin
  getPred <- function(pars, xx) pars$fmin + 
    ( (pars$f0 - pars$fmin) * pars$FractionFast) * exp( - pars$kfast * time) +
    ( (pars$f0 - pars$fmin) * (1 - pars$FractionFast) ) * exp(- pars$kslow * time)
  opar <- par()
  op<-par(mfrow=c(2,1))
  plot(all$time, all$fluorescence, xlab = "time / s", ylab = "fluorescence", main = paste(prot, conc, GEF_conc, exp, collapse = ""))
  plot(initial$time, initial$substrate_conc, xlim = c(0, 10000), xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, exp, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, kfast = 0.01, kslow = 0.001, FractionFast = 100, fmin = conc/2)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial$substrate_conc, xx = time, control = nls.lm.control(nprint = 1))
  lines(initial$time, getPred(as.list(coef(nls.out)), time), col = 2, lwd = 1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, substrate_conc = predicted.line, conc = conc)
  f0 <- coef(nls.out)[1]
  k <- coef(nls.out)[2]
  kslow <- coef(nls.out)[3]
  FractionFast <- coef(nls.out)[4]
  Plateau <- coef(nls.out)[5]
  v0 = (f0 * k * exp(k * 0)) #/ (GEF_conc*0.001)
  return(data.frame(prot, conc, f0, k, kslow, v0, FractionFast, Plateau, GEF_conc, exp))
  par(opar)
}

fit_exp_plus_linear_decay_and_plot <- function (all, initial, prot, conc, GEF_conc, exp) {
  time <- initial$time
  getPred <- function(pars, xx)  pars$plat + pars$FractionExp * ( (pars$f0 - pars$plat) * exp(- pars$k * time)) +
                              pars$f0 + (1 - pars$FractionExp) * (- pars$k_pb * time)
  opar <- par()
  op<-par(mfrow=c(2,1))
  plot(all$time, all$fluorescence, xlab = "time / s", ylab = "fluorescence", main = paste(prot, conc, GEF_conc, exp, collapse = ""))
  plot(initial$time, initial$substrate_conc, xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, exp, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(FractionExp = 0.7, f0 = conc, k = 0.01,  k_pb = 1)  # k_pb linear slope of photobleaching
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial$substrate_conc, xx = time, control = nls.lm.control(nprint = 1))
  lines(initial$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, substrate_conc = predicted.line, conc = conc)
  #plat <- coef(nls.out)[1]
  FractionExp <- coef(nls.out)[1]
  f0 <- coef(nls.out)[2]
  k <- coef(nls.out)[3]
  k_pb <- coef(nls.out)[4]
  #intersect <- coef(nls.out)[5] 
  v0 = (f0 * k * exp(k * 0))
  return(data.frame(prot, conc, f0, k, k_pb, FractionExp, v0, GEF_conc, exp))
  par(opar)
}
get_region_to_fit <- function (data) {
  if (data$substrate_conc[1] < 6) {
    plateau_fluorescence <- quantile(data$fluorescence, prob = .6, na.rm = T)
  } else {
    plateau_fluorescence <- quantile(data$fluorescence, prob = .1, na.rm = T)
  }
  plateau_times <- sort(data$time[ data$fluorescence < plateau_fluorescence ])
  cutoff_time <- plateau_times[1]
  relevant_data <- data[data$time < cutoff_time,]
  relevant_data <- relevant_data[order(relevant_data$time),]
  return(relevant_data)
}
get_linear_region <- function(data) {
  data <- cbind(data, "conc_diff" = max(data$substrate_conc) - data$substrate_conc )
  max_diff <- max(data$conc_diff)
  if (data$substrate_conc[1] < 6) {
    reaction_percentages <- seq(0.15, 0.3, 0.05)
  } else {
    reaction_percentages <- seq(0.99, 1, 0.05)
  }
  points_check <- TRUE
  for (i in reaction_percentages) {
    if (points_check) {
      linear_times <- sort(data$time[data$conc_diff < (i * max_diff)])
      if ( length (linear_times) >= 200) {
        points_check <- FALSE
      }
    }
  }
  cutoff_time <- linear_times[length(linear_times)]
  linear_data <- data[data$time < cutoff_time,]
  linear_data <- linear_data[order(linear_data$time),]
  return(linear_data)
}
fit_linear_rate <- function (data, prot, conc, GEF_conc, exp) {
  #opar <- par()
  #op<-par(mfrow=c(2,1))
  relevant.subset <- get_region_to_fit(data)
  linear.subset <- get_linear_region(relevant.subset)
  time <- unique(relevant.subset$time)
  main = paste(prot, exp, conc, GEF_conc, sep = " ")
  #y.axis.lim <- c(floor(base_fluorescence/10000)*10000, max(data$fluorescence, na.rm = T))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, main = main)
  points(relevant.subset$time, relevant.subset$fluorescence, pch = 20, col = "green")
  points(linear.subset$time, linear.subset$fluorescence, pch = 20, col = "red")
  lm.out <- lm(fluorescence ~ time, data = linear.subset)
  abline(lm(fluorescence ~ time, data = linear.subset))
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- - K
  return(data.frame(prot, conc, "rate constant" = K, "v0" = v0, GEF_conc, exp ))
  #par(opar)
}

fit_MM <- function (initial_vs_conc, prot, GEF_conc) {
  #test.conc <- seq(min(initial_vs_conc$conc), max(initial_vs_conc$conc), 0.25)
  test.conc <- unique(initial_vs_conc$conc)
  conc <- initial_vs_conc$conc
  getFit <- function(pars, xx) (conc * pars$Vmax) / (pars$Km + conc)
  residFun <- function(p, observed, xx) observed - getFit(p,xx)
  Km.estimate <- 1
  Vm.estimate <- max(initial_vs_conc$v0)
  parStart <- list(Vmax = Vm.estimate, Km = Km.estimate)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial_vs_conc$v0, xx = conc, control = nls.lm.control(nprint = 1))
  predicted.line <- getFit(as.list(coef(nls.out)), test.conc)
  predicted.line.data.frame <- data.frame(conc = test.conc, v0 = predicted.line, prot = prot)
  Vmax <- round(coef(nls.out)[1], 3)
  Km <- round(coef(nls.out)[2], 3)
  kcat = Vmax / (GEF_conc)
  standard_errors_of_parameters <- round(summary(nls.out)$coefficients[, 2], 3)
  kcat_error <- round(as.numeric(standard_errors_of_parameters[[1]]), 3) / (GEF_conc)
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
  filename = paste0("GEF_assay/", Sys.Date(), "_", prot, outsufix)
  pdf(file = filename, height = 10)
  substr.concentrations <- unique(data$conc[data$sample == prot])
  for (j in seq_along(substr.concentrations))  {
    conc = substr.concentrations[j]
    GEF_concs <- unique(data$GEF_conc[data$sample == prot & data$conc == conc])
    for (k in seq_along(GEF_concs)) {
      GEF_conc = GEF_concs[k]
      experiments <- unique(data$date[data$sample == prot & data$conc == conc & data$GEF_conc == GEF_conc])
      for (l in seq_along(experiments)) {
        exp <- experiments[l]
        data.subset <- data[ data[["sample"]] == prot & data[["conc"]] == conc & data[["GEF_conc"]] == GEF_conc & data[["date"]] == exp, c("time", "fluorescence")]
        relevant.subset <- get.substrate.conc(data.subset, conc, GEF_conc)
        #relevant.subset <- get.initial.rate.subset(data.subset, conc, GEF_conc)
        #fitting.parameters <- rbind(fitting.parameters, fit_exp_decay_and_plot(data.subset, relevant.subset, prot, conc, GEF_conc, exp))
        #fitting.parameters <- rbind(fitting.parameters, fit_two_phase_decay_and_plot(data.subset, relevant.subset, prot, conc, GEF_conc, exp))
        #fitting.parameters <- rbind(fitting.parameters, fit_exp_plus_linear_decay_and_plot(data.subset, relevant.subset, prot, conc, GEF_conc, exp))
        fitting.parameters <- rbind(fitting.parameters, fit_linear_rate(relevant.subset, prot, conc, GEF_conc, exp))
      }
    }
  }
  dev.off()
}
fitting.parameters <- cbind(fitting.parameters,paste(fitting.parameters[["GEF_conc"]], "nM", sep = " "))
names(fitting.parameters)[length(names(fitting.parameters))] <- "GEF_conc_condition"

MM.parameters <- list()
KmVm.param.df <- data.frame()
KmVm.curve.to.plot <- data.frame()
Kmkcat.df <- data.frame()
for ( i in seq_along(proteins) ) {
  mutant = proteins[i]
  GEF_concs <- unique(fitting.parameters$GEF_conc[fitting.parameters$prot == mutant])
  for (j in seq_along(GEF_concs)) {
    GEF_conc <- GEF_concs[j]
    data.subset <- fitting.parameters[fitting.parameters$prot == mutant & fitting.parameters$GEF_conc == GEF_conc,]
    MM.parameters[i] <- list (mutant = fit_MM(data.subset, mutant, GEF_conc))
    KmVm.param.df <- rbind(KmVm.param.df, MM.parameters[[i]]$parameters)
    KmVm.curve.to.plot <- rbind(KmVm.curve.to.plot, MM.parameters[[i]]$predicted_line_df)
  }
}
plot <- ggplot(data = fitting.parameters, aes(x = conc, y = v0, color = prot)) + geom_point()
plot <- plot + geom_line(data = KmVm.curve.to.plot)
#plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~Delta~'F'~s^-1))
plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~s^-1))
filename <- paste0("GEF_assay/MichaelisMenten_", Sys.Date(), ".pdf")
pdf(file = filename, width = 10)
print(plot)
dev.off()

KmVm.param.df$kcat_over_Km <- KmVm.param.df$kcat / KmVm.param.df$Km
pairwise_filename <- paste("GEF_assay/pairwise_MM_plots_", Sys.Date(), ".pdf", sep = "")
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
# 
# 
# 
# 
# 
# # fit_kobs <- function (kobs_vs_conc, prot) {
# #   test.conc <- seq(min(kobs_vs_conc), max(kobs_vs_conc), 0.25)
# #   conc <- kobs_vs_conc$conc
# #   GEF_conc <- kobs_vs_conc$GEF_conc
# #   getFit <- function (pars, xx) pars$knuc + ( (pars$K3*pars$kdiss*GEF_conc) / (1 + pars$K3*conc) )
# #   residFun <- function (p, observed, xx) observed - getFit(p,xx)
# #   parStart <- list (knuc = 0.01, K3 = 1, kdiss = 5)
# #   nls.out <- nls.lm(par = parStart, fn = residFun, observed = kobs_vs_conc$k, xx = conc, control = nls.lm.control(nprint = 1))
# #   getPred <- function(pars, xx) pars$knuc + ( (pars$K3*pars$kdiss*GEF_conc) / (1 + pars$K3*test.conc) )
# #   predicted.line <- getPred (as.list(coef(nls.out)), test.conc)
# #   predicted.line.data.frame <- data.frame(conc = test.conc, kobs = predicted.line, prot = prot)
# #   knuc <- coef(nls.out)[1]
# #   K3 <- coef(nls.out)[2]
# #   kdiss <- coef(nls.out)[3]
# #   parameters_and_predicted_line <- list("parameters" = data.frame(prot, knuc, K3, kdiss), "predicted_line_df" = predicted.line.data.frame)
# #   return(parameters_and_predicted_line)
# # }
# # kobs.parameters <- list()
# # kobs.param.df <- data.frame()
# # kobs.curve.to.plot <- data.frame()
# # kobs.df <- data.frame()
# # for (i in 1:length(proteins)) {
# #   mutant = proteins[i]
# #   data.subset <- subset(fitting.parameters, prot == mutant, select=c(conc = conc, k = k, GEF_conc = GEF_conc, exp = exp))
# #   kobs.parameters[i] <- list (mutant = fit_kobs(data.subset, mutant))
# #   kobs.param.df <- rbind(kobs.param.df, kobs.parameters[[i]]$parameters)
# #   kobs.curve.to.plot <- rbind(kobs.curve.to.plot, kobs.parameters[[i]]$predicted_line_df)
# # }
# # plot <- ggplot(data = fitting.parameters, aes(x = conc, y = k, color = prot)) + geom_point()
# # plot <- plot + geom_line(data = kobs.curve.to.plot)
# # plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('kobs / ~s^-1'))
# # filename <- paste("kobs_", Sys.Date(), ".pdf", sep = "")
# # pdf(file = filename, width = 10)
# # print(plot)
# # dev.off()
# 
#opar <- par()
filename <- paste0("GEF_assay/", Sys.Date(), "kcat_Km_GEF_exchange.pdf")
pdf(filename)
op<-par(mfrow=c(2,1))
kcat.barplot <- barplot(height = KmVm.param.df$kcat, names.arg = KmVm.param.df$prot,
                        main = "kcat of GEF exchange", ylab = bquote("kcat / s-1"),
                        xaxt = "n", cex.names = 0.75, ylim = c(0,4))
text(x = kcat.barplot, y = par("usr")[3] - 0.3, srt = 45, adj = 1, labels = KmVm.param.df$prot, xpd = TRUE)
segments(kcat.barplot, KmVm.param.df$kcat - KmVm.param.df$kcat_error, kcat.barplot,
         KmVm.param.df$kcat + KmVm.param.df$kcat_error, lwd = 1.5)
arrows(kcat.barplot, KmVm.param.df$kcat - KmVm.param.df$kcat_error, kcat.barplot,
       KmVm.param.df$kcat + KmVm.param.df$kcat_error, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
Km.barplot <- barplot(height = KmVm.param.df$Km, names.arg = KmVm.param.df$prot,
                      main = "Km of GEF exchange", ylab = bquote("Km / " *mu~M),
                      xaxt = "n", cex.names = 0.75, ylim = c(0,9))
text(x = Km.barplot, y = par("usr")[3] - 1, srt = 45, adj = 1, labels = KmVm.param.df$prot, xpd = TRUE)
segments(Km.barplot, KmVm.param.df$Km - KmVm.param.df$Km_error, Km.barplot,
          KmVm.param.df$Km + KmVm.param.df$Km_error, lwd = 1.5)
arrows(Km.barplot, KmVm.param.df$Km - KmVm.param.df$Km_error, Km.barplot,
       KmVm.param.df$Km + KmVm.param.df$Km_error, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()
op<-par(mfrow=c(1,1))
# 
# 
# 
