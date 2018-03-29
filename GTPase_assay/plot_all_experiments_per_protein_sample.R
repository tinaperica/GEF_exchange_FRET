### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
data <- read.delim("GTPase_assay/data/20170720_TP_GAP_kinetics_parsed.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc, sensor_conc, GAP_conc, gain, fluorescence)
head(data)
data$well <- do.call(paste0, data[c("row", "column")])
data$unique <- do.call(paste, c(data[, c("date", "sample", "conc", "sensor_conc", "GAP_conc", "well")], sep = " "))
data <- data[data[["conc"]] > 0,]
proteins <- as.character(unique(data$sample))
gains <- as.character(unique(data$gain))
sensor_concentrations <- unique(data$sensor_conc)
for ( i in seq_along(proteins) ) {
  plot_filename <- paste0("GTPase_assay/", proteins[i], "_GAP_assay.pdf")
  pdf(plot_filename, width = 10)
  for (j in seq_along(gains) ) {
    layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2),3,4, byrow = T))
    data_prot <- data[data$sample == proteins[i] & data$gain == gains[j],]
    unique_samples <- as.character(unique(as.character(data_prot$unique)))
    col_list <- list(unique_samples, rainbow(length(unique_samples)) )
    par(mar = c(5,4,4,1))
    plot(data_prot$time, data_prot$fluorescence, xlab = "time / s", ylab = "fluorescence", 
         type = "n") #, ylim = c(5000, 45000), xlim = c(0,10000))
    for (k in seq_along(unique_samples)) {
      temp <- data_prot[ data_prot$unique == unique_samples[k], ]
      points(temp$time, temp$fluorescence, col = col_list[[2]][k], pch = 20, cex = 0.5)
    }
    par(mar = c(5,1,4,1))
    plot(data_prot$time, data_prot$fluorescence, type = "n", axes = F, ann = F)
    legend("top", legend = col_list[[1]], title = "Ran/uM and GAP/nM", fill = col_list[[2]], bty = "n", cex = 1)
  }
  dev.off()
}
# predict.lo <- loess(fluorescence ~ time, temp, span = 0.3)
# time <- seq(0, 3000, 10)
# predicted.curve <- predict(predict.lo, time)
# plot(temp$time, temp$fluorescence)
# lines(time, predicted.curve, col = "red")
# discrets(na.omit(predict(predict.lo, time)), sym = 3)
# 
# discrets(na.omit(temp$fluorescence), sym = 3)

