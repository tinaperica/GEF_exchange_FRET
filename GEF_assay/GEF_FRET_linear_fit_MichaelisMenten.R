library(ggplot2)
library(minpack.lm)
outsufix <- paste("_fitted_curves_", Sys.Date(), ".pdf", sep = "")
tecan.data<-read.delim("GEF_assay/data/wt_check_with_sunny.txt", 
          colClasses = c('character', 'integer', 'character', 'integer', 'character', 'numeric', 'numeric', 'numeric'), head = F)
### the file with exp to discard was manually made based on inspection of all the raw data plots
#####
names(tecan.data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence)
tecan.data$well <- do.call(paste0, tecan.data[c("row", "column")])
##### Give a unique identifier to each well
#### unique is a combination of experiment date and conditions (Ran and GEF conc and well)
tecan.data<-cbind(tecan.data, paste(tecan.data$date, tecan.data$sample, tecan.data$conc, tecan.data$GEF_conc, tecan.data$well, sep = " "))
names(tecan.data)[length(names(tecan.data))] <- "unique" 

#### meaintain a manually curated file that defines wells from individual experiments to discard
data.to.discard <- read.delim("GEF_assay/data/FRET_kinetics_data_to_discard.txt", head = T) ### defined as date of exp, prot, prot conc, and GEF conc
#####
#df_args <- c(data.to.discard, sep=" ")
data.to.discard$unique <- do.call(paste, c(data.to.discard, sep=" "))
data.to.discard <- cbind(data.to.discard, paste(
  data.to.discard$exp_date, data.to.discard$protein, data.to.discard$conc, data.to.discard$GEF_conc, sep = "_"))
names(data.to.discard)[length(names(data.to.discard))] <- "unique"
#unique_identifiers_to_discard <- as.character(data.to.discard$unique)
discarded.data <- tecan.data[tecan.data$unique %in% as.character(data.to.discard$unique), ]
tecan.data <- tecan.data[! tecan.data$unique %in% as.character(data.to.discard$unique), ]
proteins<-unique(tecan.data$sample)



#### Discard the points that fit into the lowest nth percentile of fluorescence
######(this is a tunable parameter, default if 1)
##### Lowest nth percentile of data points is the linear drop due to Trp photobleaching after the reaction has finished
#### These discarded points will be plotted in black in the fitting pdf file
get_region_to_fit <- function (data) {
  plateau_fluorescence <- quantile(data$fluorescence, prob = .01, na.rm = T)
  #plateau_times <- sort(data$time[ data$fluorescence < plateau_fluorescence ])
  #cutoff_time <- plateau_times[1]
  cutoff_time <- min(data$time[ data$fluorescence < plateau_fluorescence ])
  relevant_data <- data.frame(data[data$time <= cutoff_time, ], "relevant" = "relevant")
  if (cutoff_time < max(data$time, na.rm = T)) {
    relevant_data <- rbind(relevant_data, data.frame( data[data$time > cutoff_time, ], "relevant" = "discarded") )
  }
  #relevant_data[relevant_data$relevant == "discarded",]
  return(relevant_data)
}

##### convert fluorescence to substrate (substrate is Ran:GDP, product is Ran:mant-dGDP) concentration
get.substrate.conc <- function (data, conc, GEF_conc) {
  min.f <- min(data$fluorescence[data$relevant == "relevant"], na.rm = T)
  max.f <- max(data$fluorescence[data$relevant == "relevant"], na.rm = T)
  delta.f <- max.f - min.f
  data <- cbind(data, "percent.norm.fluorescence" = 1 - (max.f - data$fluorescence)/delta.f)
  data <- cbind(data, "substrate_conc" = conc - conc*((max.f - data$fluorescence)/delta.f))
  return(data)
}


##### define a starting linear region of the data
##### Ideally that's the 10 percent of the reaction.
##### In case there are fewer than a threshold of points gradualy increase the percentage up to 30 %
get_linear_region <- function(data, conc) {
  data <- cbind(data, "conc_diff" = conc - data$substrate_conc )
  #max_diff <- max(data$conc_diff)
  max_diff <- conc
  reaction_percentages <- seq(0.10, 0.30, 0.01)
  points_check <- TRUE
  for (i in reaction_percentages) {
    if (points_check) {
      linear_times <- sort(data$time[data$conc_diff < (i * max_diff)])
      if ( length (linear_times) >= 100) {
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
  data.to.fit <- get_region_to_fit(data)
  data.converted.to.concentration <- get.substrate.conc(data.to.fit, conc, GEF_conc)
  linear.subset <- get_linear_region(data.converted.to.concentration, conc)
  simulation.time <- unique(data.converted.to.concentration$time)
  main = paste(prot, exp, conc, GEF_conc, sep = " ")
  plot(data.converted.to.concentration$time, data.converted.to.concentration$substrate_conc, 
       xlab = "time / s", ylab = "substrate concentration / uM", pch = 20, main = main)
  points(data.converted.to.concentration$time[data.converted.to.concentration$relevant == "relevant"], 
         data.converted.to.concentration$substrate_conc[data.converted.to.concentration$relevant == "relevant"], pch = 20, col = "green")
  points(linear.subset$time, linear.subset$substrate_conc, pch = 20, col = "red")
  lm.out <- lm(substrate_conc ~ time, data = linear.subset)
  abline(lm.out)
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- - K / (GEF_conc * 0.001)
  return(data.frame(prot, conc, "rate constant" = K, "v0" = v0, GEF_conc, exp ))
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
  predicted.line.data.frame <- data.frame(conc = test.conc, v0 = predicted.line, prot = prot)
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
  filename = paste0("GEF_assay/", Sys.Date(), "_", prot, outsufix)
  pdf(file = filename, height = 10)
  substr.concentrations <- unique(tecan.data$conc[tecan.data$sample == prot])
  for (j in seq_along(substr.concentrations))  {
    conc = substr.concentrations[j]
    GEF_concs <- unique(tecan.data$GEF_conc[tecan.data$sample == prot & tecan.data$conc == conc])
    for (k in seq_along(GEF_concs)) {
      GEF_conc = GEF_concs[k]
      experiments <- unique(tecan.data$date[tecan.data$sample == prot & tecan.data$conc == conc & tecan.data$GEF_conc == GEF_conc])
      for (l in seq_along(experiments)) {
        exp <- experiments[l]
        data.subset <- tecan.data[ tecan.data[["sample"]] == prot & tecan.data[["conc"]] == conc & 
                        tecan.data[["GEF_conc"]] == GEF_conc & tecan.data[["date"]] == exp, c("time", "fluorescence")]
        #relevant.subset <- get.substrate.conc(data.subset, conc, GEF_conc)
        fitting.parameters <- rbind(fitting.parameters, fit_linear_rate(data.subset, prot, conc, GEF_conc, exp))
      }
    }
  }
  dev.off()
}

fitting.parameters <- cbind(fitting.parameters,paste(fitting.parameters[["GEF_conc"]], "nM", sep = " "))
names(fitting.parameters)[length(names(fitting.parameters))] <- "GEF_conc_condition"
plot(fitting.parameters$conc, fitting.parameters$v0)

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
#filename <- paste0("GEF_assay/MichaelisMenten_", Sys.Date(), ".pdf")
#pdf(file = filename, width = 10)
print(plot)
#dev.off()
