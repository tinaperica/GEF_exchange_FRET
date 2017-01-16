library(ggplot2)
library(minpack.lm)
calibration.data <- read.delim("GTPase_assay/Sensor_calibration/20161116_Pi_titration_parsed.txt", head = F)
names(calibration.data) <- expression( experiment, time, row, column, sample, conc_string, conc, GAP_conc, fluorescence)
head(calibration.data)
times <- unique(calibration.data$time)
times/60
concentrations <- unique( calibration.data$conc)
conc_colors <- list( concentrations, rainbow(length(concentrations)) )
plot_distributions <- function (data, suffix, upylim) {
      plot(density(data$fluorescence), type = "n", 
          ylim = c(0, upylim), main = paste0("Distribution of fluorescence noise ", suffix),
          xlab = "fluorescence"
          )
      for ( i in seq_along (conc_colors[[1]]) ) {
        temp <- data[ data$conc == conc_colors[[1]][i], ]
        lines( density(temp$fluorescence), col = conc_colors[[2]][i] )
        mean.fluorescence <- mean(temp$fluorescence)
        abline(v = mean.fluorescence, col = conc_colors[[2]][i])
      }
    legend("topright", legend = conc_colors[[1]], fill = conc_colors[[2]], title = "Pi conc / uM")
}
plot_distributions(calibration.data, "over ~30 minutes", 0.002)
calibration.data.first10 <- calibration.data[ calibration.data$time < 600, ] ## data for first 10 minutes
plot_distributions(calibration.data.first10, "first 10 min", 0.002)
calibration.data.first20 <- calibration.data[ calibration.data$time < 1200, ]
plot_distributions(calibration.data.first20, "over ~20 min", 0.002)

mean.fluorescence <- with(calibration.data, aggregate(fluorescence, by=list(conc = conc), mean))
sd.fluorescence <- with(calibration.data, aggregate(fluorescence, by = list(conc = conc), sd))
fluorescence <- rbind(mean.fluorescence, data.frame(
  "conc" = mean.fluorescence$conc, "x" = mean.fluorescence$x + sd.fluorescence$x))
fluorescence <- rbind(fluorescence, data.frame(
  "conc" = mean.fluorescence$conc, "x" = mean.fluorescence$x - sd.fluorescence$x))
plot(fluorescence$conc, fluorescence$x, 
     xlab = "Pi conc / uM", ylab = "fluorescence",
     main = "0.5 uM sensor, average fluorescence over 30 min +/- 1 SD")

names(mean.fluorescence) <- c("conc", "fluorescence")
fit_data <- function (data) {
  conc <- unique(data$conc)
  #getPred <- function(pars, xx) pars$f0 + conc ( pars$b - sqrt( pars$b^2 - 4*pars$a)  )
  fmax <- max(data$fluorescence)
  f0 <- min(data$fluorescence)
  delta_f = fmax - f0
  getPred <- function(pars, xx) f0 + delta_f * (( (pars$sensor + conc + pars$Kd) - 
                    sqrt( (pars$sensor + conc + pars$Kd)^2 - 4*pars$sensor * conc ) ) / 2)
  plot(data$conc, data$fluorescence, xlab = "Pi conc / uM", ylab = "fluorescence")
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(sensor = 0.5 , Kd = 0.14)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = data$fluorescence, xx = conc, control = nls.lm.control(nprint = 1))
  lines(data$conc, getPred(as.list(coef(nls.out)), conc), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), conc)
  predicted.line.data.frame <- data.frame(conc = conc, fluorescence = predicted.line)
  sensor <- coef(nls.out)[1]
  Kd <- coef(nls.out)[2]
  conversion_factor = sensor / delta_f
  legend("bottom", legend = c(
    paste0("Kd = ", round(Kd, 2), " / uM"),
    paste0("conversion factor = ", round(conversion_factor, 8), " / uM fu-1"),
    paste0("F0 = ", round(f0, 0)),
    paste0("DeltaF = ", round(delta_f, 0))
  ), bty = "n")
  return(data.frame(f0, fmax, sensor, Kd, delta_f, conversion_factor))
}

conversion_table <- fit_data(mean.fluorescence)
conversion_table


