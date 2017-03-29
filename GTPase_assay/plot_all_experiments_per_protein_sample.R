### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
setwd("GTPase_assay")
data <- read.delim("data/GAP_kinetics.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc, GAP_conc, fluorescence)
head(data)
data$well <- do.call(paste0, data[c("row", "column")])
data$unique <- do.call(paste, c(data[, c("date", "sample", "conc", "GAP_conc", "well")], sep = " "))
data <- data[data[["conc"]] > 0,]
proteins <- as.character(unique(data$sample))
for ( i in seq_along(proteins) ) {
  plot_filename <- paste0(proteins[i], "_GAP_assay.pdf")
  pdf(plot_filename, width = 10)
  layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2),3,4, byrow = T))
  data_prot <- data[data$sample == proteins[i],]
  samples <- unique(as.character(data_prot$unique))
  col_list <- list(samples, rainbow(length(samples)) )
  par(mar = c(5,4,4,1))
  plot(data_prot$time, data_prot$fluorescence, xlab = "time / s", ylab = "fluorescence", 
       type = "n", ylim = c(5000, 45000), xlim = c(0,15000))
  for (i in seq_along(samples)) {
    temp <- data_prot[ data_prot$unique == samples[i], ]
    points(temp$time, temp$fluorescence, col = col_list[[2]][i], pch = 20, cex = 0.5)
  }
  par(mar = c(5,1,4,1))
  plot(data_prot$time, data_prot$fluorescence, type = "n", axes = F, ann = F)
  legend("top", legend = col_list[[1]], title = "Ran/uM and GAP/nM", fill = col_list[[2]], bty = "n", cex = 1)
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

