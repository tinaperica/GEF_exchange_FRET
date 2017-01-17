library(minpack.lm)
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
max_fluorescence_per_sample <- with(data, aggregate(fluorescence, by = list(sample = sample), max))
data_0.5_uM <- data[data$conc == 0.5,]
max_fluorescence_at_0.5uM <- with(data_0.5_uM, aggregate(fluorescence, by = list(sample = sample), max))
fluorescence_ratio <- merge(max_fluorescence_per_sample, max_fluorescence_at_0.5uM, by = "sample")
fluorescence_ratio <- cbind(fluorescence_ratio, data.frame("loading_efficiency" = round(fluorescence_ratio$x.y/fluorescence_ratio$x.x, 2)))
loading_efficiency <- fluorescence_ratio[,c(1,4)]
data <- merge(data, loading_efficiency, by = "sample")
data <- cbind(data, data.frame("corrected_conc" = data$conc * data$loading_efficiency))


fit_data <- function (data) {
  time <- unique(data$time)
  f0 = min(data$fluorescence, na.rm = T)
  fmax = max(data$fluorescence, na.rm = T)
  getPred <- function(pars, xx) f0 + (fmax - f0) * ( 1 - exp(- pars$K * time))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20)
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(K = 0.001)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = data$fluorescence, xx = time, control = nls.lm.control(nprint = 1))
  lines(data$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, fluorescence = predicted.line)
  K <- coef(nls.out)[1]
  v0 = (f0 * K * exp(K * 0))
  converted_rate = conversion_factor * v0
  return(data.frame("rate constant" = K, "initial_rate in fu/s" = v0, "converted_rate in uM/s" = converted_rate ))
}
temp <- data[data$sample == "PE1:GTP" & data$conc == 1 & data$date == "20170110" & data$time < 800,]
temp <- temp[order(temp$time),]
fit_data(temp)


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

