### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
data <- read.delim("GTPase_assay/data/20170110_GAP_kinetics_parsed.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc_string, conc, GAP_conc, fluorescence)
head(data)
data$prot_conc <- do.call(paste, c(data[, c("sample", "conc_string", "GAP_conc")], sep = " "))
samples <- unique(as.character(data$prot_conc))
col_list <- list(samples, rainbow(length(samples)) )
data <- data[data[["conc"]] > 0,]
plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", type = "n", ylim = c(10000, 60000))
for (i in seq_along(samples)) {
  temp <- data[ data$prot_conc == samples[i], ]
  points(temp$time, temp$fluorescence, col = col_list[[2]][i], pch = 20, cex = 0.5)
}
legend("top", legend = col_list[[1]], title = "Ran conc / uM and GAP conc / nM", fill = col_list[[2]], bty = "n", cex = 0.6)

# predict.lo <- loess(fluorescence ~ time, temp, span = 0.3)
# time <- seq(0, 3000, 10)
# predicted.curve <- predict(predict.lo, time)
# plot(temp$time, temp$fluorescence)
# lines(time, predicted.curve, col = "red")
# discrets(na.omit(predict(predict.lo, time)), sym = 3)
# 
# discrets(na.omit(temp$fluorescence), sym = 3)
