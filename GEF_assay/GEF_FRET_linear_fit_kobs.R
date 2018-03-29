library(ggplot2)
library(minpack.lm)
outsufix <- paste("_fitted_curves_kobs_", Sys.Date(), ".pdf", sep = "")
tecan.data<-read.delim(inputfilename, colClasses = c('character', 'integer', 'character', 'integer', 
                                                     'character', 'numeric', 'numeric', 'numeric', 'integer'), head = F)
### the file with exp to discard was manually made based on inspection of all the raw data plots
#####
names(tecan.data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence, cutoff_time)
tecan.data$well <- do.call(paste0, tecan.data[c("row", "column")])
##### Give a unique identifier to each well
#### unique is a combination of experiment date and conditions (Ran and GEF conc and well)
tecan.data<-cbind(tecan.data, paste(tecan.data$date, tecan.data$sample, tecan.data$conc, tecan.data$GEF_conc, tecan.data$well, sep = " "))
names(tecan.data)[length(names(tecan.data))] <- "unique" 

#### meaintain a manually curated file that defines wells from individual experiments to discard
data.to.discard <- read.delim("GEF_assay/data/FRET_kinetics_data_points_to_discard.txt", head = T) ### defined as date of exp, prot, prot conc, and GEF conc
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
get_region_to_fit <- function (total) {
  #plateau_fluorescence <- quantile(data$fluorescence, prob = .01, na.rm = T)
  #cutoff_time <- min(data$time[ data$fluorescence < plateau_fluorescence ])
  relevant_data <- data.frame(total[total$time <= total$cutoff_time, ], "relevant" = "relevant")
  relevant_data <- rbind(relevant_data, data.frame( total[total$time > total$cutoff_time, ], "relevant" = "discarded") )
  return(relevant_data)
}

get.substrate.conc <- function (total, conc, GEF_conc) {
  data <-  get_region_to_fit(total)
  min.fluor <- quantile(data$fluorescence[data$relevant == "relevant"], 0.05)
  max.fluor <- quantile(data$fluorescence[data$relevant == "relevant"], 0.9999)
  data <- cbind(data, "percent.norm.fluorescence" = (max.fluor - data$fluorescence)/(max.fluor - min.fluor))
  data <- cbind(data, "substrate_conc" = conc - conc*((max.fluor - data$fluorescence)/(max.fluor - min.fluor)))
  return(data)
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
  cutoff_time <- linear_times[length(linear_times)]
  linear_data <- data[data$time < cutoff_time,]
  linear_data <- linear_data[order(linear_data$time),]
  return(linear_data)
}
fit_linear_rate <- function (data.to.fit, prot, conc, GEF_conc, exp) {
  #data.to.fit <- get_region_to_fit(data)
  linear.subset <- get_linear_region(data.to.fit)
  main = paste(prot, exp, conc, GEF_conc, sep = " ")
  plot(data.to.fit$time, data.to.fit$substrate_conc, 
       xlab = "time / s", ylab = "substrate concentration", pch = 20, main = main)
  points(data.to.fit$time[data.to.fit$relevant == "relevant"], 
         data.to.fit$substrate_conc[data.to.fit$relevant == "relevant"], pch = 20, col = "green")
  points(linear.subset$time, linear.subset$substrate_conc, pch = 20, col = "red")
  lm.out <- lm(substrate_conc ~ time, data = linear.subset)
  abline(lm.out)
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- - K / (0.001*GEF_conc)
  return(data.frame(prot, conc, "rate constant" = K, "v0" = v0, GEF_conc, exp ))
}


fitting.parameters <- data.frame()
for (i in seq_along(proteins)) {
  prot <- proteins[i]
  filename = paste0("GEF_assay/", Sys.Date(), "_", prot, outsufix)
  pdf(file = filename, height = 10)
  substr.concentrations <- sort(as.numeric(unique(tecan.data$conc[tecan.data$sample == prot])))
  for (j in seq_along(substr.concentrations))  {
    conc = substr.concentrations[j]
    GEF_concs <- unique(tecan.data$GEF_conc[tecan.data$sample == prot & tecan.data$conc == conc])
    for (k in seq_along(GEF_concs)) {
      GEF_conc = GEF_concs[k]
      if (GEF_conc > 0) {
        experiments <- unique(tecan.data$date[tecan.data$sample == prot & tecan.data$conc == conc & tecan.data$GEF_conc == GEF_conc])
        for (l in seq_along(experiments)) {
          exp <- experiments[l]
          data.subset <- tecan.data[ tecan.data[["sample"]] == prot & tecan.data[["conc"]] == conc & 
                                   tecan.data[["GEF_conc"]] == GEF_conc & tecan.data[["date"]] == exp, c("time", "fluorescence", "cutoff_time")]
          relevant.subset <- get.substrate.conc(data.subset, conc, GEF_conc)
          fitting.parameters <- rbind(fitting.parameters, fit_linear_rate(relevant.subset, prot, conc, GEF_conc, exp))
        }
      }
    }
  }
  dev.off()
}
plot <- ggplot(data = fitting.parameters, aes(x = conc, y = v0, color = exp)) + geom_point()
print(plot)
