### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
data <- read.delim("GTPase_assay/data/GAP_kinetics.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc_string, conc, GAP_conc, fluorescence)
head(data)
samples <- as.character( unique ( data$sample ))
#samples <- c("sensor", "PE1:GTP_sensor")
col_list <- list(samples, rainbow(length(samples)) )
plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", type = "n", ylim = c(0, 5000))
for (i in seq_along( samples )) {
  temp <- data[ data$sample == samples[i], ]
  temp <- cbind( temp, "relative" = (temp$fluorescence - min(temp$fluorescence)) )
  points(temp$time, temp$fluorescence, col = col_list[[2]][i])
}
legend("top", legend = col_list[[1]], fill = col_list[[2]], bty = "n")
