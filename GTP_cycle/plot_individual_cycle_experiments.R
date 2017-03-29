library(ggplot2)

### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )
setwd("../GTP_cycle/")
exp_date <- "20170209"
data <- read.delim(paste0("data/", exp_date, "_TP_Cycle_kinetics_parsed.txt"), head = F)
names(data) <- expression( date, time, row, column, protein, prot_conc, GAP_conc, GEF_conc, MOG1_conc, fluorescence)
head(data)
data$sample <- do.call(paste, c(data[, c("protein", "prot_conc", "GAP_conc", "GEF_conc", "MOG1_conc")], sep = " "))
data <- data[data[["prot_conc"]] > 0,]
#data <- data[data$sample %in% samples,]
proteins <- as.character(unique(data$protein))
for ( i in seq_along(proteins) ) {
  plot_filename <- paste0(exp_date, "_", proteins[i], "_cycle_assay.pdf")
  pdf(plot_filename, width = 12)
  data_prot <- data[data$protein == proteins[i],]
  samples <- unique(as.character(data_prot$sample))
  col_list <- list(samples, rainbow(length(samples)) )
  samples1 <- samples[1:12]
  samples2 <- samples[13:24]
  #samples3 <- samples[10:12]
  plot <- ggplot(data = data_prot, aes(x = time, y = fluorescence, color = sample)) + geom_point(size = 1) #+ stat_smooth(span = 0.9)
  print(plot)
  plot <- ggplot(data = data_prot[data_prot$sample %in% samples1,], aes(x = time, y = fluorescence, color = sample)) + geom_point(size = 1) #+ stat_smooth(span = 0.9)
  print(plot)
  plot <- ggplot(data = data_prot[data_prot$sample %in% samples2,], aes(x = time, y = fluorescence, color = sample)) + geom_point(size = 1) #+ stat_smooth(span = 0.9)
  print(plot)
  #plot <- ggplot(data = data_prot[data_prot$sample %in% samples3,], aes(x = time, y = fluorescence, color = sample)) + geom_point(size = 0.1) + stat_smooth(span = 0.9)
  #print(plot)
  
  plot(data_prot$time, data_prot$fluorescence, xlab = "time / s", ylab = "fluorescence", type = "n", ylim = c(8000, 30000))
  for (i in seq_along(samples)) {
    temp <- data_prot[ data_prot$sample == samples[i], ]
    points(temp$time, temp$fluorescence, col = col_list[[2]][i], pch = 20, cex = 1)
  }
  legend("top", legend = col_list[[1]], title = "Ran conc/uM, GAP/uM, GEF/uM, MOG1/uM", fill = col_list[[2]], bty = "n", cex = 0.6)
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

