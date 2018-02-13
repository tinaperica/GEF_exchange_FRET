library(tidyverse)
library(minpack.lm)
##### here you define all the input and output file names
wd = "GEF_assay/" #### define working directory if needed
inputfilename <- paste0("GEF_assay/biotek_test_exp/data/good_data_parsed.txt")
#data.to.discard <- read.delim(paste0(wd, "data/FRET_kinetics_data_points_to_discard.txt"), head = T) ### defined as date of exp, prot, prot conc, and GEF conc
fit_curves_outsufix <- paste0("_fitted_curves_", ".pdf")
MM_outfilename <- paste0(wd, "output_plots/MichaelisMenten_", Sys.Date(), ".pdf")

biotek.data <- read_delim(inputfilename, delim = "\t", col_names = T)

proteins <- unique(biotek.data$sample)

( data.to.fit <- biotek.data %>% 
  group_by(condition) %>%
  mutate(relevant = ifelse(Time < cutoff_time, "relevant", "discarded")) %>%
  filter(relevant == "relevant") %>%
  mutate(max.fluorescence = quantile(fluorescence, prob = 0.99, na.rm = T), 
         min.fluorescence = quantile(fluorescence, prob = 0.01, na.rm = T)) %>%
  mutate(substrate_conc = conc - conc*((max.fluorescence - fluorescence)/(max.fluorescence - min.fluorescence)))
)
conditions <- unique(data.to.fit$condition)

fit_linear_decay <- function (data) {
  simulation.time <- data$Time
  condition <- data$condition[1]
  conc <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  linear.subset <- data %>% filter(substrate_conc > 0.75*conc)
  lm.out <- lm(substrate_conc ~ Time, data = linear.subset)
  f0 <- coefficients(lm.out)[1]
  K <- coefficients(lm.out)[2]
  v0 <- - K / (0.001*GEF_conc)
  return(data.frame("v0" = v0, "intercept" = f0, "slope" = K))
}

fit_exp_decay <- function (data) {
  simulation.time <- data$Time
  getPred <- function(pars, xx) (pars$f0 - pars$c) * exp(- pars$k * simulation.time) + pars$c
  condition <- data$condition[1]
  conc <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  residFun <- function(p, observed, xx) observed - getPred(p,xx)
  parStart <- list(f0 = conc, k = 0.01, c = 0)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = data$substrate_conc, xx = simulation.time, control = nls.lm.control(nprint = 1))
  predicted.line <- getPred(as.list(coef(nls.out)), simulation.time)
  predicted.line.df <- tibble("time" = simulation.time, "conc" = predicted.line)
  k <- coef(nls.out)[2]
  f0 <- coef(nls.out)[1]
  v0 = (f0 * k * exp(k * 0)) / (GEF_conc*0.001)  #### this is not real v0, this is v0 normalized for GEF concentration
  return(tibble("v0" = v0, "predicted.conc" = predicted.line.df$conc, "predicted.time" = predicted.line.df$time))
}

fit_MM <- function (fitted.parameters) {
  test.conc <- unique(fitted.parameters$conc)
  conc <- fitted.parameters$conc
  GEF_conc <- fitted.parameters$GEF_conc[1]
  getFit <- function(pars, xx) (conc * pars$Vmax) / (pars$Km + conc)
  residFun <- function(p, observed, xx) observed - getFit(p,xx)
  Km.estimate <- 1
  Vm.estimate <- max(fitted.parameters$v0)
  parStart <- list(Vmax = Vm.estimate, Km = Km.estimate)
  nls.out <- nls.lm(par = parStart, fn = residFun, observed = fitted.parameters$v0, xx = conc, control = nls.lm.control(nprint = 1))
  predicted.line <- getFit(as.list(coef(nls.out)), conc)
  predicted.line.data.frame <- data.frame(conc = conc, v0 = predicted.line)
  Vmax <- round(coef(nls.out)[1], 3)
  Km <- round(coef(nls.out)[2], 3)
  kcat = Vmax
  standard_errors_of_parameters <- round(summary(nls.out)$coefficients[, 2], 3)
  kcat_error <- round(as.numeric(standard_errors_of_parameters[[1]]), 3) / (GEF_conc)
  Km_error <- round(as.numeric(standard_errors_of_parameters[[2]]), 3)
  summary <- summary(nls.out)
  sd_of_fit <- round(as.numeric(summary[["sigma"]]), 3)
  return(tibble("Vmax" = Vmax,  "Km" = Km, "kcat" = kcat, "kcat_error" = kcat_error, "Km_error" = Km_error,
                "sd_of_fit" = sd_of_fit, "predicted.conc" = predicted.line.data.frame$conc, "predicted.line" = predicted.line.data.frame$v0))
}

exp.fitting.parameters <-  data.to.fit %>% group_by(sample, conc, condition, GEF_conc) %>%
  #filter(conc < 10) %>%
  select(Time, substrate_conc, condition, sample, conc, GEF_conc) %>%
  do(data.frame(fit_exp_decay(.)))

plots <- list()
for (i in seq_along(conditions)) {
  fitted.data <- exp.fitting.parameters %>% filter(condition == conditions[i])
  concentrations <- unique(fitted.data$conc)
  raw.data <- data.to.fit %>% filter(condition == conditions[i] & conc %in% concentrations)
  if (nrow(raw.data) > 0) {
    plots[[i]] <- ggplot(raw.data, mapping = aes(x = Time, y = substrate_conc)) + 
      geom_point(color = "black") + geom_line(fitted.data, mapping = aes(x = predicted.time, y = predicted.conc), color = "red") + 
      ggtitle(conditions[i])
  }
}

pdf("GEF_assay/exp_fitted_data.pdf")
print(plots)
dev.off()

# linear.fitting.parameters <-  data.to.fit %>% group_by(sample, conc, condition, GEF_conc) %>%
#   select(Time, substrate_conc, condition, sample, conc, GEF_conc) %>%
#   filter(conc >= 12) %>%
#   do(data.frame(fit_linear_decay(.)))
# mean.linear.fitting.parameters <- linear.fitting.parameters %>% 
#   mutate("round_conc" = round(conc, 0)) %>%
#   group_by(sample, round_conc) %>%
#   summarise("mean_v0" = mean(v0))
# 
# plots <- list()
# for (i in seq_along(conditions)) {
#   fitted.data <- linear.fitting.parameters %>% filter(condition == conditions[i])
#   concentrations <- fitted.data$conc
#   raw.data <- data.to.fit %>% filter(condition == conditions[i] & conc %in% concentrations)
#   if (nrow(raw.data) > 0) {
#     plots[[i]] <- ggplot(raw.data, mapping = aes(x = Time, y = substrate_conc)) + 
#       geom_point(color = "black") + geom_abline(slope = fitted.data$slope, intercept = fitted.data$intercept, color = "red") + 
#       ggtitle(conditions[i])
#   }
# }
# 
# pdf("GEF_assay/linear_fitted_data.pdf")
# print(plots)
# dev.off()

#ggplot(linear.fitting.parameters, aes(x = conc, y = v0, color = sample)) + geom_point()

#ggplot(mean.linear.fitting.parameters, aes(x = round_conc, y = mean_v0, color = sample)) + geom_point()

#chosen.linear.fitting.parameters <- linear.fitting.parameters %>% 
 #   select(sample, conc, condition, GEF_conc, v0) %>%
  #  filter(conc >= 8)
#chosen.exp.fitting.parameters <- exp.fitting.parameters %>%
 # select(sample, conc, condition, GEF_conc, v0) %>%
  #filter(conc < 8)
  
#combined.fitting.parameters <- bind_rows(chosen.exp.fitting.parameters, chosen.linear.fitting.parameters)

#combined.fitting.parameters <- bind_rows(exp.fitting.parameters, linear.fitting.parameters)
combined.fitting.parameters <- exp.fitting.parameters

mean.combined.fitting.parameters <- combined.fitting.parameters %>% 
  mutate("round_conc" = round(conc, 1)) %>%
  group_by(sample, round_conc) %>%
  summarise("mean_v0" = mean(v0), "plus_err" = mean_v0 + var(v0), "minus_err" = mean_v0 - var(v0))

MM.parameters <-  combined.fitting.parameters %>% group_by(sample) %>%
  select(sample, conc, v0, GEF_conc) %>%
  do(data.frame(fit_MM(.)))

ggplot(combined.fitting.parameters, aes(x = conc, y = v0, color = sample)) + 
  geom_point() + geom_line(MM.parameters, mapping = aes(x = predicted.conc, y = predicted.line, color = sample))

ggplot(mean.combined.fitting.parameters, aes(x = round_conc, y = mean_v0, color = sample)) + 
        geom_point() + geom_smooth(MM.parameters, mapping = aes(x = predicted.conc, y = predicted.line, color = sample))




