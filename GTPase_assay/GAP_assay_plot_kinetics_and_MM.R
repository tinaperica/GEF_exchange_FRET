### plot and analyse GTPase assay with phosphate sensor data
options( stringsAsFactors = F )

conversion_factor = 2.760382e-05 # uM / fluorescence
base_fluorescence = 12471.8

data <- read.delim("GTPase_assay/data/GAP_kinetics.txt", head = F)
names(data) <- expression( date, time, row, column, sample, conc_string, conc, GAP_conc, fluorescence)
head(data)
#### convert fluorescence into Pi concentration
data <- cbind(data, data.frame("delta_f" = data$fluorescence - base_fluorescence))

### correct substrate concentration based on loading
proteins <- as.character(unique(data$sample))
data$prot_conc <- do.call(paste, c(data[, c("sample", "conc_string")], sep = " "))
samples <- unique(as.character(data$prot_conc))
col_list <- list(samples, rainbow(length(samples)) )
data <- data[data[["conc"]] > 0,]
plot(data$time, data$Pi_conc, xlab = "time / s", ylab = "fluorescence", type = "n")
for (i in seq_along(samples)) {
    temp <- data[ data$prot_conc == samples[i], ]
    #temp <- cbind( temp, "relative" = (temp$fluorescence - min(temp$fluorescence)) )
    points(temp$time, temp$Pi_conc, col = col_list[[2]][i], pch = 20, cex = 0.5)
  }
legend("top", legend = col_list[[1]], title = "Ran conc / uM", fill = col_list[[2]], bty = "n", cex = 0.6)

