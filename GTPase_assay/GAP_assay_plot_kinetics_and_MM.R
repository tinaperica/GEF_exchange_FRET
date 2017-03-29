library(minpack.lm)
library(ggplot2)

### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
setwd("GTPase_assay")
conversion_factor = 0.0001064242 # uM / fluorescence
base_fluorescence = 9726.956
fluorescence.at.2uM.Pi <- 27173.39 - base_fluorescence
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")

data <- read.delim("data/GAP_kinetics.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc, GAP_conc, fluorescence)
head(data)
data$well <- do.call(paste0, data[c("row", "column")])

data.to.discard <- read.delim("data/GAP_assay_data_to_discard.txt", head = T) ### defined as date of exp, prot, prot conc, and GEF conc
#####
data<-cbind(data, paste(data$date, data$sample, data$conc, data$GAP_conc, data$well, sep = " "))
names(data)[length(names(data))] <- "unique"  #### unique is a combination of experiment date and conditions (Ran and GAP conc and well)

df_args <- c(data.to.discard, sep=" ")
data.to.discard$unique <- do.call(paste, df_args)


unique_identifiers_to_discard <- as.character(data.to.discard$unique)
data.discarded <- data[data$unique %in% unique_identifiers_to_discard,]
data <- data[! data$unique %in% unique_identifiers_to_discard,]

### correct substrate concentration based on loading
#max_fluorescence_per_sample <- with(data, aggregate(fluorescence, by = list(sample = sample), max))
#expected_fluorescence_at_2uM <- data.frame("sample" = as.character(unique(data$sample)), "x" = fluorescence.at.2uM.Pi) 
#data_at_2_uM <- data[data$conc == 2,]
#max_fluorescence_at_2_uM <- with(data_at_2_uM, aggregate(fluorescence, by = list(sample = sample), max))
#loaded_conc_from_max_fluorescence_at_2_uM <- cbind(max_fluorescence_at_2_uM, "conc_loaded" = max_fluorescence_at_2_uM$x * conversion_factor)
#loading_efficiency <- cbind(loaded_conc_from_max_fluorescence_at_2_uM, "loading_efficiency" = loaded_conc_from_max_fluorescence_at_2_uM$conc_loaded/2)
#data <- merge(data, loading_efficiency, by = "sample")
#data <- cbind(data, data.frame("corrected_conc" = data$conc * data$loading_efficiency))
# ##
data <- cbind(data, data.frame("corrected_conc" = data$conc))

proteins <- as.character(unique(data$sample))
proteins <- proteins[proteins != "buffer"]
get_region_to_fit <- function (data) {
 plateau_fluorescence <- quantile(data$fluorescence, prob = .95, na.rm = T)
 plateau_times <- sort(data$time[ data$fluorescence > plateau_fluorescence ])
 cutoff_time <- plateau_times[1] + 50
 relevant_data <- data[data$time < cutoff_time,]
 relevant_data <- relevant_data[order(relevant_data$time),]
 return(relevant_data)
}
get_linear_region <- function(data) {
  relevant_data <- get_region_to_fit(data)
  relevant_data <- cbind(relevant_data, "fluor_diff" = relevant_data$fluorescence - min(relevant_data$fluorescence))
  linear_times <- sort(relevant_data$time[relevant_data$fluor_diff < (0.7 * max(relevant_data$fluor_diff))])
  cutoff_time <- linear_times[-1]
  relevant_data <- data[data$time < cutoff_time,]
  relevant_data <- relevant_data[order(relevant_data$time),]
  return(relevant_data)
}
fit_one_phase_association <- function (data) {
  opar <- par()
  op<-par(mfrow=c(2,1))
  relevant.subset <- get_region_to_fit(data)
  time <- unique(relevant.subset$time)
  #f0 = min(relevant.subset$fluorescence, na.rm = T)
  #fmax = max(relevant.subset$fluorescence, na.rm = T)
  getPred <- function(pars, xx) pars$f0 + (pars$fmax - pars$f0) * ( 1 - exp(- pars$K * time))
  main = paste(data$sample[1], data$date[1], data$corrected_conc[1], data$GAP_conc[1], sep = " ")
  y.axis.lim <- c(floor(base_fluorescence/10000)*10000, max(data$fluorescence, na.rm = T))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, main = main, ylim = y.axis.lim)
  plot(relevant.subset$time, relevant.subset$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, ylim = y.axis.lim)
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(K = 0.001, 
                   f0 = min(relevant.subset$fluorescence, na.rm = T),
                   fmax = max(relevant.subset$fluorescence, na.rm = T))
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = relevant.subset$fluorescence, xx = time, control = nls.lm.control(nprint = 1))
  lines(relevant.subset$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, fluorescence = predicted.line)
  K <- coef(nls.out)[1]
  f0 <- coef(nls.out)[2]
  fmax <- coef(nls.out)[3]
  v0 = ( (fmax - f0) * K * exp(- K * 0) )
  return(data.frame("rate constant" = K, "initial_rate" = v0 ))
  par(opar)
}


fit_linear_rate <- function (data) {
  #opar <- par()
  #op<-par(mfrow=c(2,1))
  relevant.subset <- get_linear_region(data)
  time <- unique(relevant.subset$time)
  #f0 = min(relevant.subset$fluorescence, na.rm = T)
  #fmax = max(relevant.subset$fluorescence, na.rm = T)
  main = paste(data$sample[1], data$date[1], data$corrected_conc[1], data$GAP_conc[1], sep = " ")
  y.axis.lim <- c(floor(base_fluorescence/10000)*10000, max(data$fluorescence, na.rm = T))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, main = main, ylim = y.axis.lim)
  points(relevant.subset$time, relevant.subset$fluorescence, pch = 20, col = "red")
  lm.out <- lm(fluorescence ~ time, data = relevant.subset)
  abline(lm(fluorescence ~ time, data = relevant.subset))
  #predicted.line.data.frame <- abline(lm(fluorescence ~ time, data = relevant.subset))
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- K
  return(data.frame("rate constant" = K, "initial_rate" = v0 ))
  #par(opar)
}

# SpanFast=(Plateau-Y0)*PercentFast*.01
# SpanSlow=(Plateau-Y0)*(100-PercentFast)*.01
# Y=Y0+ SpanFast*(1-exp(-KFast*X)) + SpanSlow*(1-exp(-KSlow*X))


fit_MichaelisMenten <- function (initial_vs_conc, prot) {
  test.conc <- seq(min(initial_vs_conc$conc), max(initial_vs_conc$conc), 0.25)
  conc <- sort(initial_vs_conc$conc)
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

fitting_parameters <- data.frame()
for (p in seq_along(proteins)) {
  prot <- proteins[p]
  #filename = paste0(Sys.Date(), "_", substr(prot, 1,4), outsufix)
  filename = paste0(Sys.Date(), "_", prot, outsufix)
  pdf(file = filename, height = 10)
  corrected_concentrations <- unique(data$corrected_conc[data$sample == prot])
  for (cc in seq_along(corrected_concentrations)) {
    conc <- corrected_concentrations[cc]
    GAP_concentrations <- unique(data$GAP_conc[data$sample == prot & data$corrected_conc == conc])
    for (g in seq_along(GAP_concentrations)) {
      GAP_conc = GAP_concentrations[g]
      experiments <- unique(data$date[data$sample == prot & data$corrected_conc == conc & data$GAP_conc == GAP_conc])
      for (e in seq_along(experiments)) {
        exp <- experiments[e]
        wells <- unique(data$well[data$sample == prot & data$corrected_conc == conc & data$GAP_conc == GAP_conc & data$date == exp])
        for (w in seq_along(wells)) {
          well <- wells[w]
          data.subset <- data[data$sample == prot & data$corrected_conc == conc & data$GAP_conc == GAP_conc & data$date == exp & data$well == well,]
          data.subset <- data.subset[order(data.subset$time),]
          #rates <- fit_linear_rate(data.subset)
          rates <- fit_one_phase_association(data.subset)
          rates <- cbind(rates, data.frame(prot, conc, GAP_conc, exp, well))
          fitting_parameters <- rbind(fitting_parameters, rates)
        }
      }
    }
  }
  dev.off()
}
fitting_parameters <- cbind(fitting_parameters, data.frame("converted_rate" = fitting_parameters$initial_rate * conversion_factor))
fitting_parameters <- cbind(fitting_parameters, data.frame(
  "v0" = ifelse(fitting_parameters$GAP_conc > 0,
              fitting_parameters$converted_rate / (fitting_parameters$GAP_conc/1000),
              fitting_parameters$converted_rate)
  ))
fitting_parameters <- fitting_parameters[fitting_parameters$GAP_conc != 0,]
MM.parameters <- list()
KmVm.param.df <- data.frame()
KmVm.curve.to.plot <- data.frame()
Kmkcat.df <- data.frame()
for (i in seq_along(proteins)) {
  mutant = proteins[i]
  data.subset <- fitting_parameters[fitting_parameters$prot == mutant, ]
  data.subset <- data.subset[, c("conc", "v0", "GAP_conc")]
  data.subset <- data.subset[order(data.subset$conc),]
  MM.parameters[i] <- list (mutant = fit_MichaelisMenten(data.subset, mutant))
  KmVm.param.df <- rbind(KmVm.param.df, MM.parameters[[i]]$parameters)
  KmVm.curve.to.plot <- rbind(KmVm.curve.to.plot, MM.parameters[[i]]$predicted_line_df)
}
plot <- ggplot(data = fitting_parameters, aes(x = conc, y = v0, color = prot)) + geom_point()
plot <- plot + geom_line(data = KmVm.curve.to.plot)
#plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote("Initial exchange rate " ~v[0]~ " /  [GAP conc]" ~s^-1))

plot <- plot + xlab(expression('Ran:GTP concentration / ' * mu ~ M)) + ylab(expression('Initial exchange rate' ~ v[0] ~ '/[' ~ mu ~ M ~ 'GAP conc]' ~ s^-1))

filename <- paste0("MichaelisMenten_", Sys.Date(), ".pdf")
pdf(file = filename, width = 10)
print(plot)
dev.off()


KmVm.param.df <- KmVm.param.df[KmVm.param.df$prot != "PE3:GTP",]

Km.barplot <- barplot(height = KmVm.param.df$Km, names.arg = KmVm.param.df$prot, 
                      main = "Km of GAP mediated GTP hydrolysis", ylab = bquote("Km / " *mu~M),
                      xaxt = "n", cex.names = 0.75, ylim = c(0,10))
text(x = Km.barplot, y = par("usr")[3] - 1, srt = 45, adj = 1, labels = KmVm.param.df$prot, xpd = TRUE)
segments(Km.barplot, KmVm.param.df$Km - KmVm.param.df$Km_error, Km.barplot,
         KmVm.param.df$Km + KmVm.param.df$Km_error, lwd = 1.5)
arrows(Km.barplot, KmVm.param.df$Km - KmVm.param.df$Km_error, Km.barplot,
       KmVm.param.df$Km + KmVm.param.df$Km_error, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)



kcat.barplot <- barplot(height = KmVm.param.df$kcat, names.arg = KmVm.param.df$prot, 
                        main = "kcat of GAP mediated GTP hydrolysis", ylab = bquote("kcat / s-1"),
                        xaxt = "n", cex.names = 0.6, ylim = c(0,0.35))
text(x = kcat.barplot, y = par("usr")[3] - 0.3, srt = 45, adj = 1, labels = KmVm.param.df$prot, xpd = TRUE)
segments(kcat.barplot, KmVm.param.df$kcat - KmVm.param.df$kcat_error, kcat.barplot,
         KmVm.param.df$kcat + KmVm.param.df$kcat_error, lwd = 1.5)
arrows(kcat.barplot, KmVm.param.df$kcat - KmVm.param.df$kcat_error, kcat.barplot,
       KmVm.param.df$kcat + KmVm.param.df$kcat_error, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)




