library(ggplot2)

conversion_factor <- 3.684e-5
base_fluorescence <- 18158
data <- read.delim("data/20170224_TP_Cycle_kinetics_parsed.txt", head = F)
names(data) <- expression( date, time, row, column, protein, prot_conc, GAP_conc, GEF_conc, MOG1_conc, fluorescence)

get_region_to_fit <- function (data) {
  plateau_fluorescence <- quantile(data$fluorescence, prob = .95, na.rm = T)
  plateau_times <- sort(data$time[ data$fluorescence > plateau_fluorescence ])
  cutoff_time <- plateau_times[1] + 500
  relevant_data <- data[data$time < cutoff_time,]
  relevant_data <- relevant_data[order(relevant_data$time),]
  return(relevant_data)
}
fit_one_phase_association <- function (data) {
  opar <- par()
  op<-par(mfrow=c(2,1))
  relevant.subset <- get_region_to_fit(data)
  time <- unique(relevant.subset$time)
  f0 = min(relevant.subset$fluorescence, na.rm = T)
  fmax = max(relevant.subset$fluorescence, na.rm = T)
  getPred <- function(pars, xx) f0 + (fmax - f0) * ( 1 - exp(- pars$K * time))
  main = paste(data$sample[1], data$date[1], data$corrected_conc[1], 
        data$GAP_conc[1], data$GEF_conc[1], data$MOG1_conc[1], sep = " ")
  y.axis.lim <- c(floor(base_fluorescence/10000)*10000, max(data$fluorescence, na.rm = T))
  plot(data$time, data$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, main = main, ylim = y.axis.lim)
  plot(relevant.subset$time, relevant.subset$fluorescence, xlab = "time / s", ylab = "fluorescence", pch = 20, ylim = y.axis.lim)
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(K = 0.001)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = relevant.subset$fluorescence, xx = time, control = nls.lm.control(nprint = 1))
  lines(relevant.subset$time, getPred(as.list(coef(nls.out)), time), col=2, lwd=1)
  predicted.line <- getPred(as.list(coef(nls.out)), time)
  predicted.line.data.frame <- data.frame(time = time, fluorescence = predicted.line)
  K <- coef(nls.out)[1]
  v0 = ( (fmax - f0) * K * exp(- K * 0) )
  return(data.frame("rate constant" = K, "initial_rate" = v0 ))
  par(opar)
}
fitting_parameters <- data.frame()
filename <- paste0(Sys.Date(), "_cycle_fitted.pdf")
pdf(filename)
proteins <- unique(data$protein)
GAP_concs <- unique(data$GAP_conc)
GEF_concs <- unique(data$GEF_conc)
MOG1_concs <- unique(data$MOG1_conc) 
prot_concs <- unique(data$prot_conc)
data_to_test_plot <- list()
for (protein in proteins) {
  for (gap in seq_along( GAP_concs ) ) {
    for (gef in seq_along( GEF_concs ) ) {
        data_to_test_plot[[paste(GAP_concs[gap], GEF_concs[gef], sep = "_")]] <- data[
          data$protein == protein & data$GAP_conc == GAP_concs[gap] & data$GEF_conc == GEF_concs[gef],]
        for (mog1 in seq_along( MOG1_concs ) ) {
          for (prot in seq_along( prot_concs )) {
            data_sub <- data[data$GAP_conc == GAP_concs[gap] & 
                         data$GEF_conc == GEF_concs[gef] &
                          data$prot_conc == prot_concs[prot] &
                         data$MOG1_conc == MOG1_concs[mog1], ]
            if ( length (data_sub[,1]) > 0 ) { 
              rates <- fit_one_phase_association(data_sub)
              rates <- cbind( rates, data.frame("protein" = protein, "protein_conc" = protein_conc,
                      "GAP_conc" = GAP_concs[gap],
                      "GEF_conc" = GEF_concs[gef], "MOG1_conc" = MOG1_concs[mog1]),
                      "max" = max_fluorescence <- max(data_sub$fluorescence, na.rm = T) )
              fitting_parameters <- rbind(fitting_parameters, rates)
          }
        }
      }
    }
  }
}
dev.off()
pdf(paste(Sys.Date(), "GTP_cycle_vs_MOG1_conc_raw_data.pdf", sep = "_"), width = 10)
plot <- ggplot(data = data_to_test_plot[["0.02_0.01"]][data_to_test_plot[["0.02_0.01"]]$time < 7500,], aes(x = time, y = fluorescence, color = MOG1_conc)) + coord_cartesian(ylim=c(7500,28000)) + geom_point()
plot <- plot + ggtitle("2 uM Gsp1, 20 nM GAP and 10 nM GEF")
print(plot)
plot <- ggplot(data = data_to_test_plot[["0.02_0.1"]][data_to_test_plot[["0.02_0.1"]]$time < 7500,], aes(x = time, y = fluorescence, color = MOG1_conc)) + coord_cartesian(ylim=c(7500,28000)) + geom_point()
plot <- plot + ggtitle("2 uM Gsp1, 20 nM GAP and 100 nM GEF")
print(plot)
plot <- ggplot(data = data_to_test_plot[["0.1_0.01"]][data_to_test_plot[["0.1_0.01"]]$time < 7500,], aes(x = time, y = fluorescence, color = MOG1_conc)) + coord_cartesian(ylim=c(7500,28000)) + geom_point()
plot <- plot + ggtitle("2 uM Gsp1, 100 nM GAP and 10 nM GEF")
print(plot)
dev.off()
fitting_parameters <- cbind(fitting_parameters, data.frame("converted_rate" = fitting_parameters$initial_rate * conversion_factor))
fitting_parameters <- cbind(fitting_parameters, data.frame(
  "v0" = ifelse(fitting_parameters$GAP_conc > 0,
                fitting_parameters$converted_rate / fitting_parameters$GAP_conc,
                fitting_parameters$converted_rate)
))
fitting_parameters$gap_gef <- do.call(paste, c(fitting_parameters[, c("GAP_conc", "GEF_conc")], sep = "_"))

pdf(paste(Sys.Date(), "rate_vs_MOG1_conc.pdf", sep = "_"), width = 10) 
plot <- ggplot (data = fitting_parameters, aes(x = MOG1_conc, y = v0, color = gap_gef)) + geom_point() + stat_smooth(level = 0.9)
plot <- plot + xlab(bquote("MOG1 concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~'/ [GAP conc]' ~s^-1))
print(plot)
plot <- ggplot (data = fitting_parameters, aes(x = MOG1_conc, y = converted_rate, color = gap_gef)) + geom_point() + stat_smooth(level = 0.9)
plot <- plot + xlab(bquote("MOG1 concentration / " *mu~M)) + ylab(bquote('Initial exchange rate '~v[0]~ ~s^-1))
print(plot)
plot <- ggplot (data = fitting_parameters, aes(x = MOG1_conc, y = max, col = gap_gef)) + geom_point() + stat_smooth(level = 0.9)
plot <- plot + xlab(bquote("MOG1 concentration / " *mu~M)) + ylab(bquote("Max fluorescence"))
print(plot)
dev.off()

write.table(file = "MOG1_GAP_cycle_rates.txt", fitting_parameters[,c(3:6, 8)], quote = F, row.names = F, sep = "\t")
