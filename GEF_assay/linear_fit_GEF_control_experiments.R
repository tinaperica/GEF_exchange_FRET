library(ggplot2)

### plot and fit GEF control experiments
options( stringsAsFactors = F )
outsufix <- paste("_control_exp_fitted_", Sys.Date(), ".pdf", sep = "")

data <- read.delim("GEF_assay/data/20170424_TP_FRET_control_kinetics_parsed.txt", head = F)
head(data)
names(data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence, nucleotide)
data$well <- do.call(paste0, data[c("row", "column")])
data$condition <- do.call(paste, data[c("GEF_conc", "nucleotide")])

control_data <- data[data$conc == 3,]

proteins <- as.character(unique(control_data$sample))
get_region_to_fit <- function (data) {
  plateau_fluorescence <- quantile(data$fluorescence, prob = .7, na.rm = T)
  plateau_times <- sort(data$time[ data$fluorescence < plateau_fluorescence ])
  cutoff_time <- plateau_times[1]
  relevant_data <- data[data$time < cutoff_time,]
  relevant_data <- relevant_data[order(relevant_data$time),]
  return(relevant_data)
}

get.substrate.conc <- function (total, conc) {
  min.total.fluor <- quantile(total$fluorescence, 0.001)
  max.total.fluor <- quantile(total$fluorescence, 0.999)
  total <- cbind(total, "percent.norm.fluorescence" = (max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor))
  total <- cbind(total, "substrate_conc" = conc - conc*((max.total.fluor - total$fluorescence)/(max.total.fluor - min.total.fluor)))
  return(total)
}
get_region_to_fit <- function (data) {
  plateau_fluorescence <- quantile(data$fluorescence, prob = .5, na.rm = T)
  plateau_times <- sort(data$time[ data$fluorescence < plateau_fluorescence ])
  cutoff_time <- plateau_times[1]
  relevant_data <- data[data$time < cutoff_time,]
  relevant_data <- relevant_data[order(relevant_data$time),]
  return(relevant_data)
}
get_linear_region <- function(data, min_percent, max_percent, n_points) {
  data <- cbind(data, "conc_diff" = max(data$substrate_conc) - data$substrate_conc )
  max_diff <- max(data$conc_diff)
  reaction_percentages <- seq(min_percent, max_percent, 0.05)
  points_check <- TRUE
  for (i in reaction_percentages) {
    if (points_check) {
      linear_times <- sort(data$time[data$conc_diff < (i * max_diff)])
      if ( length (linear_times) >= n_points) {
        points_check <- FALSE
      }
    }
  }
  cutoff_time <- linear_times[length(linear_times)]
  linear_data <- data[data$time < cutoff_time,]
  linear_data <- linear_data[order(linear_data$time),]
  return(linear_data)
}
fit_linear_rate <- function (data, prot, conc, GEF_conc, nucleotide, exp) {
  if (GEF_conc > 0) {
    relevant.subset <- get_region_to_fit(data)
    linear.subset <- get_linear_region(relevant.subset, 0.1, 0.5, 100)
  } else {
    linear.subset <- data
  }
  time <- unique(relevant.subset$time)
  main = paste(prot, exp, conc, GEF_conc, nucleotide, sep = " ")
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, main = main)
  points(relevant.subset$time, relevant.subset$fluorescence, pch = 20, col = "green")
  points(linear.subset$time, linear.subset$fluorescence, pch = 20, col = "red")
  lm.out <- lm(fluorescence ~ time, data = linear.subset)
  abline(lm(fluorescence ~ time, data = linear.subset))
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- - K
  return(data.frame(prot, conc, "rate constant" = K, "v0" = v0, GEF_conc, nucleotide, exp ))
}


fitted.data <- data.frame()
for (i in seq_along(proteins)) {
  prot <- proteins[i]
  filename = paste0("GEF_assay/", Sys.Date(), "_", prot, outsufix)
  pdf(file = filename, width = 10)
  substr.concentrations <- unique(control_data$conc[control_data$sample == prot])
  per_protein_data <- control_data[control_data$sample == prot,]
  per_protein_plot <- ggplot(data = per_protein_data, aes(x = time, y = fluorescence, color = condition)) + geom_point()
  print(per_protein_plot)
  for (j in seq_along(substr.concentrations))  {
    conc = substr.concentrations[j]
    GEF_concs <- unique(control_data$GEF_conc[control_data$sample == prot & control_data$conc == conc])
    for (k in seq_along(GEF_concs)) {
      GEF_conc = GEF_concs[k]
      nucleotides <- unique(control_data$nucleotide[control_data$sample == prot & control_data$conc == conc & control_data$GEF_conc == GEF_conc])
      for (l in seq_along(nucleotides)) {
        nucleotide <- nucleotides[l]
        experiments <- unique(control_data$date[control_data$sample == prot & control_data$conc == conc 
                            & control_data$GEF_conc == GEF_conc & control_data$nucleotide == nucleotide])
        for (m in seq_along(experiments)) {
          exp <- experiments[m]
          data.subset <- control_data[ control_data[["sample"]] == prot & control_data[["conc"]] == conc & control_data[["GEF_conc"]] == GEF_conc & 
                          control_data[["nucleotide"]] == nucleotide & control_data[["date"]] == exp, c("time", "fluorescence")]
          relevant.subset <- get.substrate.conc(data.subset, conc)
          fitted.data <- rbind(fitted.data, fit_linear_rate(relevant.subset, prot, conc, GEF_conc, nucleotide, exp))
        }
      }
    }
  }
  dev.off()
}

fitted.data$condition <- do.call(paste, fitted.data[c("GEF_conc", "nucleotide")])

p<-ggplot(fitted.data, aes(x=prot, y=v0, color = condition)) +
  geom_bar(stat="identity", position=position_dodge(), fill="white")
print(p)
