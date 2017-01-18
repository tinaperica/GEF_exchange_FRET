library(minpack.lm)
library(ggplot2)

### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )

conversion_factor = 2.760382e-05 # uM / fluorescence
base_fluorescence = 12471.8
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")

data <- read.delim("GTPase_assay/data/GAP_kinetics.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc_string, conc, GAP_conc, fluorescence)
head(data)
data$well <- do.call(paste0, data[c("row", "column")])
data <- data[data$GAP_conc <= 30,]
### correct substrate concentration based on loading
max_fluorescence_per_sample <- with(data, aggregate(fluorescence, by = list(sample = sample), max))
data_0.5_uM <- data[data$conc == 0.5,]
max_fluorescence_at_0.5uM <- with(data_0.5_uM, aggregate(fluorescence, by = list(sample = sample), max))
fluorescence_ratio <- merge(max_fluorescence_per_sample, max_fluorescence_at_0.5uM, by = "sample")
fluorescence_ratio <- cbind(fluorescence_ratio, data.frame("loading_efficiency" = round(fluorescence_ratio$x.y/fluorescence_ratio$x.x, 2)))
loading_efficiency <- fluorescence_ratio[,c(1,4)]
data <- merge(data, loading_efficiency, by = "sample")
data <- cbind(data, data.frame("corrected_conc" = data$conc * data$loading_efficiency))
##
proteins <- as.character(unique(data$sample))

get_region_to_fit <- function (data) {
 plateau_fluorescence <- quantile(data$fluorescence, prob = .95, na.rm = T)
 plateau_times <- sort(data$time[ data$fluorescence > plateau_fluorescence ])
 cutoff_time <- plateau_times[1] + 200
 relevant_data <- data[data$time < cutoff_time,]
 relevant_data <- relevant_data[order(relevant_data$time),]
 return(relevant_data)
}
fit_one_phase_association <- function (data) {
  opar <- par()
  op<-par(mfrow=c(2,1))
  relevant.subset <- get_region_to_fit(data)
  time <- unique(relevant.subset$time)
  f0 = min(relevant.subset$fluorescence, na.rm = T)
  fmax = max(relevant.subset$fluorescence, na.rm = T)
  getPred <- function(pars, xx) f0 + (fmax - f0) * ( 1 - exp(- pars$K * time))
  main = paste(data$sample[1], data$date[1], data$corrected_conc[1], data$GAP_conc[1], sep = " ")
  y.axis.lim <- c(floor(base_fluorescence/10000)*10000, max(data$fluorescence, na.rm = T))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, main = main, ylim = y.axis.lim)
  plot(relevant.subset$time, relevant.subset$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, ylim = y.axis.lim)
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(K = 0.001)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = relevant.subset$fluorescence, xx = time, control = nls.lm.control(nprint = 1))
  lines(relevant.subset$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, fluorescence = predicted.line)
  K <- coef(nls.out)[1]
  v0 = (fmax - f0) * K * exp(- K * 0)
  return(data.frame("rate constant" = K, "initial_rate" = v0 ))
  par(opar)
}
fit_MichaelisMenten <- function (initial_vs_conc, prot) {
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

fitting_parameters <- data.frame()
for (p in seq_along(proteins)) {
  prot <- proteins[p]
  filename = paste0(Sys.Date(), "_", substr(prot, 1,4), outsufix)
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
fitting_parameters <- cbind(fitting_parameters, data.frame("v0" = fitting_parameters$converted_rate / (fitting_parameters$GAP_conc/1000)))
plot <- ggplot(data = fitting_parameters, aes(x = conc, y = v0, color = prot)) + geom_point()
print(plot)
