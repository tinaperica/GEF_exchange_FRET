library(ggplot2)
library(minpack.lm)
options( stringsAsFactors = F )
data <- read.delim("GTPase_assay/Sensor_calibration/20170717_TP_sensor_calibration_parsed.txt", head = F)
names(data) <- expression( experiment, time, row, column, sample, conc, sensor_conc, fluorescence)
head(data)
data <- data[order(data$conc, data$time),]
times <- unique(data$time)
sensor_concentrations <- unique(data$sensor_conc)
conversion_table <- data.frame()

plot_distributions <- function (data, suffix, upylim) {
  data <- data[complete.cases(data),]
  concentrations <- unique(data$conc)
  conc_colors <- list( concentrations, rainbow(length(concentrations)) )
  plot(density(data$fluorescence), type = "n", 
       ylim = c(0, upylim), main = paste0("Distribution of fluorescence noise ", suffix),
       xlab = "fluorescence"
  )
  for ( i in seq_along (conc_colors[[1]]) ) {
    temp <- data[ data$conc == conc_colors[[1]][i], ]
    if ( length(temp[,1]) > 3 ) {
      lines( density(temp$fluorescence), col = conc_colors[[2]][i] )
      mean.fluorescence <- mean(temp$fluorescence)
      abline(v = mean.fluorescence, col = conc_colors[[2]][i])
    }
  }
  legend("topright", legend = conc_colors[[1]], fill = conc_colors[[2]], title = "Pi conc / uM")
}

fit_data <- function (data, data.with.sd) {
  #conc <- unique(data$conc)
  #sensor_conc <- 4
  #### equation to fit
  #getPred <- function(pars, xx) pars$f0 + pars$delta_f * (( (pars$sensor_conc + conc + pars$Kd) - 
  #sqrt( (pars$sensor_conc + conc + pars$Kd)^2 - 4*pars$sensor_conc * conc ) ) / (2*pars$sensor_conc) )
  #plot(data$conc, data$fluorescence, xlab = "Pi conc / uM", ylab = "fluorescence")
  #residFun <- function(p, observed, xx) observed - getPred(p,xx)
  #starting parameters
  #parStart <- list(f0 = 7000, delta_f = 45000, Kd = 0.1, sensor_conc = 10)
  #nls.out <- nls.lm(par = parStart, fn = residFun, observed = data$fluorescence, xx = conc, control = nls.lm.control(nprint = 1))
  #lines(data$conc, getPred(as.list(coef(nls.out)), conc), col=2, lwd=1)
  #predicted.line <- getPred(as.list(coef(nls.out)), conc)
  #predicted.line.data.frame <- data.frame(conc = conc, fluorescence = predicted.line)
  #fluorescence.at.2uM.Pi <- predicted.line.data.frame$fluorescence[predicted.line.data.frame$conc == 2]
  #f0 <- coef(nls.out)[1]
  #delta_f <- coef(nls.out)[2]
  #Kd <- coef(nls.out)[3]
  #sensor_conc <- round(coef(nls.out)[4], 2)
  #conversion_factor = sensor_conc / delta_f
  plot(data.with.sd$conc, data.with.sd$fluorescence, 
       xlab = "Pi conc / uM", ylab = "fluorescence",
       main = paste0(sensor_conc, " uM sensor, average fluorescence +/- 1 SD"))
  linear_range <- data[data$conc < 0.6*sensor_conc, ]
  linear_fit <- lm(fluorescence ~ conc, data = linear_range)
  abline(linear_fit, col = "green")
  linear_slope <- coefficients(linear_fit)[2]
  linear_intercept <- coefficients(linear_fit)[1]
  legend("bottom", legend = c(
    #paste0("Kd = ", round(Kd, 2), " / uM"),
    #paste0("conversion factor = ", round(conversion_factor, 8), " / uM fu-1"),
    #paste0("F0 = ", round(f0, 0)),
    #paste0("DeltaF = ", round(delta_f, 0)),
    paste0("sensor conc. / uM = ", sensor_conc),
    paste0("linear slope = ", linear_slope),
    paste0("intercept = ", linear_intercept)
  ), bty = "n")
  #return(data.frame(f0, sensor_conc, Kd, delta_f, conversion_factor, linear_slope, linear_intercept))
  return(data.frame(sensor_conc, linear_slope, linear_intercept))
  
}
mean_fluorescence_all <- data.frame()
for ( sc in seq_along(sensor_concentrations) ) {
  sensor_conc <- sensor_concentrations[sc]
  calibration.data <- data[data$sensor_conc == sensor_conc, ]
  signal_distribution_filename <- paste0("GTPase_assay/Sensor_calibration/", 
        Sys.Date(), "_", sensor_conc, "_distribution_of_fluorescence_signal.pdf")
  pdf(signal_distribution_filename)
  plot_distributions(calibration.data, "over entire experiment", 0.002)
  dev.off()
  mean.fluorescence <- with(calibration.data, aggregate(fluorescence, by=list(conc = conc), mean))
  sd.fluorescence <- with(calibration.data, aggregate(fluorescence, by = list(conc = conc), sd))
  names(mean.fluorescence) <- c("conc", "fluorescence")
  mean_fluorescence_all <- rbind(mean_fluorescence_all, data.frame("sensor" = as.character(sensor_conc), mean.fluorescence))
  fluorescence.plus.sd <- rbind(mean.fluorescence, data.frame(
    "conc" = mean.fluorescence$conc, "fluorescence" = mean.fluorescence$fluorescence + sd.fluorescence$x))
  fluorescence.minus.sd <- rbind(mean.fluorescence, data.frame(
    "conc" = mean.fluorescence$conc, "fluorescence" = mean.fluorescence$fluorescence - sd.fluorescence$x))
  fluorescence.with.sd <- rbind(fluorescence.plus.sd, fluorescence.minus.sd)
  linear_fit_filename <- paste0("GTPase_assay/Sensor_calibration/",
            Sys.Date(), "_linear_sensor_calibration_", sensor_conc, ".pdf")
  pdf(linear_fit_filename)
  conversion_table <- rbind(conversion_table, fit_data(mean.fluorescence, fluorescence.with.sd))
  dev.off()
  
}

conversion_table
mean_fluorescence_all
pdf(file = paste("GTPase_assay/Sensor_calibration/", Sys.Date(), "mean_fluorescence_all_sensor_conc.pdf"))
plot <- ggplot(data = mean_fluorescence_all, aes(x = conc, y = fluorescence, color = sensor)) + geom_point() + geom_abline()
print(plot)
dev.off()
