### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
setwd("GTPase_assay")
exp_date <- "20170327"
data <- read.delim(paste0("data/", exp_date, "_TP_GAP_kinetics_parsed.txt"), head = F)
names(data) <- expression( date, time, row, column, sample, conc, GAP_conc, fluorescence)
head(data)
data$prot_conc <- do.call(paste, c(data[, c("date", "sample", "conc", "GAP_conc", "row", "column")], sep = " "))
data <- data[data[["conc"]] > 0,]
proteins <- as.character(unique(data$sample))
for ( i in seq_along(proteins) ) {
  plot_filename <- paste0(exp_date, "_", proteins[i], "_GAP_assay.pdf")
  pdf(plot_filename)
  data_prot <- data[data$sample == proteins[i],]
  samples <- unique(as.character(data_prot$prot_conc))
  col_list <- list(samples, rainbow(length(samples)) )
  plot(data_prot$time, data_prot$fluorescence, xlab = "time / s", ylab = "fluorescence", type = "n", ylim = c(5000, 50000))
  for (i in seq_along(samples)) {
    temp <- data_prot[ data_prot$prot_conc == samples[i], ]
    points(temp$time, temp$fluorescence, col = col_list[[2]][i], pch = 20, cex = 0.5)
  }
  legend("top", legend = col_list[[1]], title = "Ran conc / uM and GAP conc / nM", fill = col_list[[2]], bty = "n", cex = 0.6)
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

