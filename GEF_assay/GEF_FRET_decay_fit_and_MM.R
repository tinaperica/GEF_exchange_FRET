library(ggplot2)
library(minpack.lm)
##### here you define all the input and output file names
wd = "GEF_assay/" #### define working directory if needed
inputfilename <- paste0("GEF_assay/biotek_test_exp/data/good_data_parsed.txt")
data.to.discard <- read.delim(paste0(wd, "data/FRET_kinetics_data_points_to_discard.txt"), head = T) ### defined as date of exp, prot, prot conc, and GEF conc
fit_curves_outsufix <- paste0("_fitted_curves_", ".pdf")
MM_outfilename <- paste0(wd, "output_plots/MichaelisMenten_", Sys.Date(), ".pdf")

#tecan.data<-read.delim(inputfilename, colClasses = c('character', 'integer', 'character', 'integer', 
 #                         'character', 'numeric', 'numeric', 'numeric', 'integer'), head = F)
tecan.data <- read.delim(inputfilename, head = T)
### the file with exp to discard was manually made based on inspection of all the raw data plots
#####
#names(tecan.data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence, cutoff_time)
tecan.data$well <- do.call(paste0, tecan.data[c("row", "column")])
##### Give a unique identifier to each well
#### unique is a combination of experiment date and conditions (Ran and GEF conc and well)
tecan.data<-cbind(tecan.data, paste(tecan.data$date, tecan.data$sample, tecan.data$conc, tecan.data$GEF_conc, tecan.data$well, sep = " "))
names(tecan.data)[length(names(tecan.data))] <- "unique" 
tecan.data <- tecan.data[order(tecan.data$conc), ]
#### meaintain a manually curated file that defines wells from individual experiments to discard
##### named data_to_discard.txt
#df_args <- c(data.to.discard, sep=" ")
data.to.discard$unique <- do.call(paste, c(data.to.discard, sep=" "))
data.to.discard <- cbind(data.to.discard, paste(
  data.to.discard$exp_date, data.to.discard$protein, data.to.discard$conc, data.to.discard$GEF_conc, sep = "_"))
names(data.to.discard)[length(names(data.to.discard))] <- "unique"
#unique_identifiers_to_discard <- as.character(data.to.discard$unique)
discarded.data <- tecan.data[tecan.data$unique %in% as.character(data.to.discard$unique), ]
tecan.data <- tecan.data[! tecan.data$unique %in% as.character(data.to.discard$unique), ]
proteins<-unique(tecan.data$sample)


##### convert fluorescence to substrate (substrate is Ran:GDP, product is Ran:mant-dGDP) concentration
get.substrate.conc <- function (data, conc, GEF_conc) {
  relevant_data <- data.frame(data[data$time < data$cutoff_time, ], "relevant" = "relevant")
  discarded.data <- data[data$time > data$cutoff_time, ]
  if (length( discarded.data[,1]) > 0) {
    relevant_data <- rbind(relevant_data, data.frame(data[data$time >= data$cutoff_time, ], "relevant" = "discarded"))
  }
  min.f <- min(relevant_data$fluorescence[relevant_data$relevant == "relevant"], na.rm = T)
  max.f <- max(relevant_data$fluorescence[relevant_data$relevant == "relevant"], na.rm = T)
  delta.f <- max.f - min.f
  relevant_data <- cbind(relevant_data, "percent.norm.fluorescence" = 1 - (max.f - relevant_data$fluorescence)/delta.f)
  relevant_data <- cbind(relevant_data, "substrate_conc" = conc - conc*((max.f - relevant_data$fluorescence)/delta.f))
  return(relevant_data)
}

fit_exp_decay_and_plot <- function (data, prot, conc, GEF_conc, exp, date) {
  data.to.fit <- data[data$relevant == "relevant", ]
  simulation.time <- data.to.fit$time
  getPred <- function(pars, xx) (pars$f0 - pars$c) * exp(- pars$k * simulation.time) + pars$c
  opar <- par()
  op<-par(mfrow=c(2,1))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", main = paste(prot, conc, GEF_conc, exp, collapse = ""))
  plot(data.to.fit$time, data.to.fit$substrate_conc, xlab = "time / s", ylab = "substrate concentration", main = paste(prot, conc, GEF_conc, exp, collapse = " "))
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, k = 0.01, c = 0)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = data.to.fit$substrate_conc, xx = simulation.time, control = nls.lm.control(nprint = 1))
  lines(data.to.fit$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=3)
  predicted.line <- getPred(as.list(coef(nls.out)), simulation.time)
  predicted.line.data.frame <- data.frame(time = simulation.time, substrate_conc = predicted.line, conc = conc)
  k <- coef(nls.out)[2]
  f0 <- coef(nls.out)[1]
  v0 = (f0 * k * exp(k * 0)) / (GEF_conc*0.001)  #### this is not real v0, this is v0 normalized for GEF concentration
  return(data.frame(prot, conc, f0, k, v0, GEF_conc, exp, date))
  par(opar)
}

fit_MM <- function (initial_vs_conc, prot, GEF_conc) {
  test.conc <- unique(initial_vs_conc$conc)
  conc <- initial_vs_conc$conc
  getFit <- function(pars, xx) (conc * pars$Vmax) / (pars$Km + conc)
  residFun <- function(p, observed, xx) observed - getFit(p,xx)
  Km.estimate <- 1
  Vm.estimate <- max(initial_vs_conc$v0)
  parStart <- list(Vmax = Vm.estimate, Km = Km.estimate)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial_vs_conc$v0, xx = conc, control = nls.lm.control(nprint = 1))
  predicted.line <- getFit(as.list(coef(nls.out)), test.conc)
  predicted.line.data.frame <- data.frame(conc = conc, v0 = predicted.line, prot = prot)
  Vmax <- round(coef(nls.out)[1], 3)
  Km <- round(coef(nls.out)[2], 3)
  kcat = Vmax
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
  filename = paste0(wd, "output_plots/", Sys.Date(), "_", prot, fit_curves_outsufix)
  pdf(file = filename, height = 10)
  substr.concentrations <- unique(tecan.data$conc[tecan.data$sample == prot])
  for (j in seq_along(substr.concentrations))  {
    conc = substr.concentrations[j]
    GEF_concs <- unique(tecan.data$GEF_conc[tecan.data$sample == prot & tecan.data$conc == conc])
    for (k in seq_along(GEF_concs)) {
      GEF_conc = GEF_concs[k]
      if (GEF_conc > 0) {
        experiments <- unique(tecan.data$unique[tecan.data$sample == prot & tecan.data$conc == conc & tecan.data$GEF_conc == GEF_conc])
        for (l in seq_along(experiments)) {
          exp <- experiments[l]
          date <- tecan.data[["date"]][ tecan.data[["unique"]] == exp][1]
          data.subset <- tecan.data[ tecan.data[["sample"]] == prot & tecan.data[["conc"]] == conc & 
                                     tecan.data[["GEF_conc"]] == GEF_conc & tecan.data[["unique"]] == exp, c("Time", "fluorescence", "cutoff_time"),]
          data.subset <- get.substrate.conc(data.subset, conc, GEF_conc)
          fitting.parameters <- rbind(fitting.parameters, fit_exp_decay_and_plot(data.subset, prot, conc, GEF_conc, exp, date))
        }
      }
    }
  }
  dev.off()
}

#fitting.parameters <- cbind(fitting.parameters,paste(fitting.parameters[["GEF_conc"]], "nM", sep = " "))
#names(fitting.parameters)[length(names(fitting.parameters))] <- "GEF_conc_condition"
plot <- ggplot(fitting.parameters[fitting.parameters$v0 < 20,], aes(x = conc, y = v0, color = date)) + geom_point()#, xlab = "substrate concentration", ylab = "v0 / GEF conc in uM")
print(plot)

write.table(file = paste0(wd, "fitting_parameters.txt"), fitting.parameters, quote = F, row.names = F, sep = "\t")

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
plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~s^-1))
pdf(file = MM_outfilename, width = 10)
print(plot)
dev.off()
write.table(file = paste0(wd, "MichaelisMenten_parameters.txt"), KmVm.param.df, quote = F, row.names = F, sep = "\t")
