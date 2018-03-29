library(ggplot2)
##### here you define all the input and output file names
wd = "GEF_assay/" #### define working directory if needed
inputfilename <- paste0(wd, "data/GEF_FRET_kinetics.txt")
#inputfilename <- paste0(wd, "data/20170329_TP_FRET_kinetics_parsed.txt")

data.to.discard <- read.delim(paste0(wd, "data/FRET_kinetics_data_points_to_discard.txt"), head = T) ### defined as date of exp, prot, prot conc, and GEF conc
outsufix <- paste0("_fitted_curves_", ".pdf")
MM_outfilename <- paste0(wd, "output_plots/MichaelisMenten_", Sys.time(), ".pdf")

tecan.data<-read.delim(inputfilename, colClasses = c('character', 'integer', 'character', 'integer', 
                                                     'character', 'numeric', 'numeric', 'numeric', 'integer'), head = F)
### the file with exp to discard was manually made based on inspection of all the raw data plots
#####
names(tecan.data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence, cutoff_time)
tecan.data <- tecan.data[tecan.data$GEF_conc == 5,]
max.time <- max( tecan.data[["time"]], na.rm = T) #### use to have all the plots on the same time scale
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


get.substrate.conc <- function (data, conc, GEF_conc) {
  relevant_data <- data.frame(data[data$time < data$cutoff_time, ], "relevant" = "relevant")
  if ( length( data[data$time > data$cutoff_time, ][,1] ) > 0 ) {
    relevant_data <- rbind(relevant_data, data.frame(data[data$time >= data$cutoff_time, ], "relevant" = "discarded"))
  }
  min.f <- quantile(relevant_data$fluorescence[relevant_data$relevant == "relevant"], 0.1)
  max.f <- max(relevant_data$fluorescence[relevant_data$relevant == "relevant"], 0.9)
  delta.f <- max.f - min.f
  relevant_data <- cbind(relevant_data, "percent.norm.fluorescence" = 1 - (max.f - relevant_data$fluorescence)/delta.f)
  relevant_data <- cbind(relevant_data, "substrate_conc" = conc - conc*((max.f - relevant_data$fluorescence)/delta.f))
  return(relevant_data)
}


##### define a starting linear region of the data
##### Ideally that's the 10 percent of the reaction.
##### In case there are fewer than a threshold of points gradualy increase the percentage up to 30 %
get_linear_region <- function(data) {
  max.f <- max(data$fluorescence)
  data <- cbind(data, "diff" = max.f - data$fluorescence)
  max_diff <- max(data$fluorescence) - min(data$fluorescence)
  reaction_percentages <- seq(0.10, 0.5, 0.01)
  points_check <- TRUE
  for (i in reaction_percentages) {
    if (points_check) {
      linear_times <- sort(data$time[data$diff < (i * max_diff)])
      if ( length (linear_times) >= 50) {
        points_check <- FALSE
      }
    }
  }
  linear_cutoff_time <- linear_times[length(linear_times)]
  linear_data <- data[data$time < linear_cutoff_time,]
  linear_data <- linear_data[order(linear_data$time),]
  return(linear_data)
}
fit_linear_rate <- function (data, prot, conc, GEF_conc, exp, date) {
  linear.subset <- get_linear_region(data[ data[["relevant"]] == "relevant", ])
  main = paste(prot, exp, conc, GEF_conc, sep = " ")
  xlim <- c(0, max.time)
  plot(data$time, data$substrate_conc, 
       xlab = "time / s", ylab = "substrate concentration", pch = 20, xlim = xlim, main = main)
  plot(data$time[data$relevant == "relevant"], data$substrate_conc[data$relevant == "relevant"], 
       xlab = "time / s", ylab = "substrate concentration", pch = 20, xlim = xlim)
  points(linear.subset$time, linear.subset$substrate_conc, pch = 20, col = "red")
  lm.out <- lm(substrate_conc ~ time, data = linear.subset)
  abline(lm.out)
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- - K / (0.001*GEF_conc)  ### v0 is divided by the enzyme concentration in uM
  return(data.frame(prot, conc, "rate constant" = K, "v0" = v0, GEF_conc, exp, date ))
}

fit_MM <- function (initial_vs_conc, prot) {
  conc <- initial_vs_conc$conc
  getFit <- function(pars, xx) (conc * pars$Vmax) / (pars$Km + conc)
  residFun <- function(p, observed, xx) observed - getFit(p,xx)
  Km.estimate <- 1
  Vm.estimate <- max(initial_vs_conc$v0)
  parStart <- list(Vmax = Vm.estimate, Km = Km.estimate)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = initial_vs_conc$v0, xx = conc, control = nls.lm.control(nprint = 1))
  predicted.line <- getFit(as.list(coef(nls.out)), conc)
  predicted.line.data.frame <- data.frame(conc = conc, v0 = predicted.line, prot = prot)
  Vmax <- round(coef(nls.out)[1], 3)
  Km <- round(coef(nls.out)[2], 3)
  kcat = Vmax
  standard_errors_of_parameters <- round(summary(nls.out)$coefficients[, 2], 3)
  kcat_error <- round(as.numeric(standard_errors_of_parameters[[1]]), 3) / (GEF_conc)
  Km_error <- round(as.numeric(standard_errors_of_parameters[[2]]), 3)
  summary <- summary(nls.out)
  sd_of_fit <- round(as.numeric(summary[["sigma"]]), 3)
  predicted.line.data.frame <- cbind(predicted.line.data.frame, 
                data.frame("max" = predicted.line.data.frame$v0 + sd_of_fit, "min" = predicted.line.data.frame$v0 - sd_of_fit))
  parameters_and_predicted_line <- list("parameters" = data.frame(prot, Vmax, Km, kcat, kcat_error, Km_error, sd_of_fit), 
                                        "predicted_line_df" = predicted.line.data.frame)
  return(parameters_and_predicted_line)
}

fitting.parameters <- data.frame()
for (i in seq_along(proteins)) {
  prot <- proteins[i]
  filename = paste0("GEF_assay/output_plots/", Sys.time(), "_", prot, outsufix)
  pdf(file = filename, height = 10)
  opar <- par()
  op <- par(mfrow = c(2, 1))
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
                                       tecan.data[["GEF_conc"]] == GEF_conc & tecan.data[["unique"]] == exp, c("time", "fluorescence", "cutoff_time")]
          relevant.subset <- get.substrate.conc(data.subset, conc, GEF_conc)
          fitting.parameters <- rbind(fitting.parameters, fit_linear_rate(relevant.subset, prot, conc, GEF_conc, exp, date))
        }
      }
    }
  }
  par(opar)
  dev.off()
}

plot <- ggplot(fitting.parameters, aes(x = conc, y = v0, color = date)) + geom_point()#, xlab = "substrate concentration", ylab = "v0 / GEF conc in uM")
print(plot)
write.table(file = paste0(wd, "fitting_parameters.txt"), fitting.parameters, quote = F, row.names = F, sep = "\t")

MM.parameters <- list()
KmVm.param.df <- data.frame()
KmVm.curve.to.plot <- data.frame()
Kmkcat.df <- data.frame()
for ( i in seq_along(proteins) ) {
  mutant = proteins[i]
  data.subset <- fitting.parameters[ fitting.parameters$prot == mutant, ]
  MM.parameters[i] <- list (mutant = fit_MM(data.subset, mutant))
  KmVm.param.df <- rbind(KmVm.param.df, MM.parameters[[i]]$parameters)
  KmVm.curve.to.plot <- rbind(KmVm.curve.to.plot, MM.parameters[[i]]$predicted_line_df)
}

plot <- ggplot(data = fitting.parameters, aes(x = conc, y = v0, color = prot)) + geom_point()
plot <- plot + geom_line(data = KmVm.curve.to.plot)
plot <- plot + geom_ribbon(data = KmVm.curve.to.plot, aes(ymin = min, ymax = max), linetype=2, alpha=0.1)
plot <- plot + xlab(bquote("Substrate (Ran:GDP) concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GEF conc]' ~s^-1))

pdf(file = MM_outfilename, width = 10)
print(plot)
dev.off()
write.table(file = paste0(wd, "MichaelisMenten_parameters.txt"), KmVm.param.df, quote = F, row.names = F, sep = "\t")
