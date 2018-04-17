# import libraries
library(tidyverse)
library(lubridate)
library(minpack.lm)



# read in data
inputfilename <- "GEF_assay/biotek_test_exp/data/good_data_parsed.txt"
biotek.data <- read_delim(inputfilename, delim = "\t", col_names = T)


run_nls <- function(data, deadtime, debug = FALSE) {
  
  c0 <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  
  #### debug mode: if debug is TRUE, print out condition being fit 
  if (debug) { print(unique(data$condition)) }
  
  #### add the (estimated) dead time
  if (deadtime > 0) {
    deadtime_timepoints <- seq(-1*(deadtime), -5, 5)
    timepoints_to_append <- data[1:length(deadtime_timepoints),]
    timepoints_to_append$Time <- deadtime_timepoints
    timepoints_to_append$observed <- NA
    data <- rbind(data, timepoints_to_append)
  }
  start <- list(span1 = max(data$observed, na.rm = T) - min(data$observed, na.rm = T),
                f0 = max(data$observed, na.rm = T),
                f_plateau = 0.7 * max(data$observed, na.rm = T),
                k = 0.001,
                k_background = 3e-5)
  out <- nlsLM(observed ~ span1 * exp(-k * Time) + (f0 - span1 - f_plateau) * exp(-k_background * Time) + f_plateau,
                      data = data, start = start, control = nls.lm.control(maxiter = 300))
  #out <- nlsLM(observed ~ (f0 - f_plateau * exp(-k_background * Time)) * exp(-k * Time) + 
   #              f_plateau * exp(-k_background * Time) + 5e5 * (k - 3e-5),
    #           data = data, start = start, control = nls.lm.control(maxiter = 300))
  
  #### save optimal parameters
  span1 <- coef(out)[1]
  f0 <- coef(out)[2]
  f_plateau <- coef(out)[3]
  k <- coef(out)[4]
  k_background <- coef(out)[5]
  span2 <- f0 - span1 - f_plateau
  
  if (k_background < 0) {
    start <- list(f0 = max(data$observed, na.rm = T),
                  k = 0.001,
                  #c = median(data$observed, na.rm = T),
                  f_plateau = 0.7 * max(data$observed, na.rm = T)
                  )
    out <- nlsLM(observed ~ (f0 - f_plateau) * exp(-k * Time) + f_plateau,
                 data = data, start = start, control = nls.lm.control(maxiter = 300))
    f0 <- coef(out)[1]
    k <- coef(out)[2]
    f_plateau <- coef(out)[3]
    k_background <- 0
    span1 <- f0 - f_plateau
    span2 <- 0
  }
  #### calculate ideal curves for overall fit, true exchange signal, and background, save
  #data$predicted <- (f0 - f_plateau * exp(-k_background * data$Time)) * exp(-k * data$Time) + f_plateau * exp(-k_background * data$Time)
  #data$exchange <- (f0 - f_plateau) * (exp(-k * data$Time)) + f_plateau
  #data$background <- f_plateau * exp(-k_background * data$Time)
  data$predicted <- span1 * exp(-k * data$Time) +  span2* exp(-k_background * data$Time) + f_plateau
  data$exchange <- span1 * (exp(-k * data$Time)) + f_plateau
  data$background <-  span2 * exp(-k_background * data$Time) + f_plateau
  data$k_background <- k_background
  
  #### save optimal parameters in the data table
  data$f0 <- f0
  data$k <- k
  data$f_plateau <- f_plateau
  data$k_background <- k_background
  data$span1 <- span1
  data$span2 <- span2
  ### calculate the initial rate
  data$vf0 = (f0 * k * exp(k * 0)) / (GEF_conc*0.001)   ### initial rate in fluorescence units
  data$v0 <- (c0 * k * exp(k * 0)) / (GEF_conc*0.001)   ## initial rate in concentration units  
  data$f_signal_change <- (f0 - f_plateau) / f0 
  return(data)
}

plot_data <- function(data, debug = FALSE) {
  
  #### Call this to get rid of the "[[1]] [[2]] [[3]]..." printing during pdf plotting. note that it clears errors, must be removed when debugging
  sink('/dev/null')  
  
  #### set date, for filenaming purposes
  today = gsub('-', '', today(tzone="US/Pacific"))
  
  #### plot data with fitted curves
  plots <- list()
  conditions = unique(data$condition)
  for (i in seq_along(conditions)) {
    data_to_plot <- data %>% filter(condition == conditions[i])
    plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) + geom_point(color = "black") +
      geom_line(data_to_plot, mapping = aes(x = Time, y = predicted), color = "red") +
      geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
      geom_line(data_to_plot, mapping = aes(x = Time, y = background), color = "blue") +
      ggtitle(conditions[i])
  }
  pdf(paste0(today, 'data.pdf'))
  print(plots)
  dev.off()
  
  #### make michaelis-menten plot of v0 vs. conc
  # plots <- list()
  # samples = unique(data$sample)
  # for (i in seq_along(samples)) {
  #   data_to_plot <- data %>% filter(sample == samples[i])
  #   plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = conc, y = v0, color = as.character(date))) + geom_point() +
  #     ggtitle(samples[i])
  # }
  # pdf(paste0(today, 'MM.pdf'), width = 10)
  # print(plots)
  # dev.off()
  
  #### plot histograms showing distribution of parameters across all fits
  if (debug) {
    pdf(paste(today, 'params.pdf', sep = '_'))
    #print(ggplot(data, aes(x = conc, y = f0, color = sample)) + geom_point())
    print(ggplot(data, aes(x = conc, y = k, color = sample)) + geom_point())
    #print(ggplot(data, aes(x = conc, y = f_plateau, color = sample)) + geom_point())
    #print(ggplot(data, aes(x = conc, y = f_signal_change, color = sample)) + geom_point())
    print(ggplot(data, aes(x = conc, y = k_background, color = sample)) + geom_point())
    print(ggplot(data, aes(x = conc, y = span1, color = sample)) + geom_point())
    print(ggplot(data, aes(x = conc, y = span2, color = sample)) + geom_point())
    #print(ggplot(data, aes(x = conc, y = v0, color = sample)) + geom_point())
    #print(ggplot(data, aes(x = conc, y = vf0, color = sample)) + geom_point())
    dev.off()
  }
  
  return(data)
}

fit_MM <- function(data) {
  start <- list(Vmax = 5, Km = 2)
  out <- nlsLM(v0 ~ (conc * Vmax) / (Km + conc),
               data = data, start = start, control = nls.lm.control(maxiter = 200))
  Vmax <- coef(out)[1]
  Km <- coef(out)[2]
  kcat <- Vmax
  data$predicted_v0 <- (data$conc * Vmax) / (Km + data$conc)
  data$kcat <- kcat
  data$Km <- Km
  return(data)
}

#### Make a copy of data to be transformed
data.exp <- biotek.data

#### Fit data assuming photobleaching decay is exponential, with observations in fluorescence units
data.exp <- data.exp %>%
  group_by(condition) %>%
  mutate(observed = fluorescence) %>%
  do(run_nls(., deadtime = 18, debug = T))
fit.parameters <- data.exp %>% 
  select(sample, date, conc, condition, f0, k, k_background, f_plateau, vf0, v0, f_signal_change) %>% 
  unique() %>%
  arrange(sample, conc)
ggplot(fit.parameters, aes(conc, v0, color = sample)) + geom_point()

data.exp <- plot_data(data.exp, debug = T)



samples <- fit.parameters %>% pull(sample) %>% unique()
#samples <- samples[c(1:4,6)]
MM.plots <- list()
for (s in seq_along(samples)) {
  temp <- fit.parameters %>% filter(sample == samples[s])
  MM.plots[[s]] <- ggplot(temp, aes(x = conc, y = v0, color = as.character(date))) + 
    geom_point() + ggtitle(samples[s])
}
pdf("MM_curve.pdf", width = 10)
print(MM.plots)
dev.off()



MM.data <- tibble(v0 = double(), conc = double(), 
                sample = character(), predicted_v0 = double(), 
                kcat = double(), Km = double())
for (samp in samples) {
  temp <- fit.parameters %>%
    ungroup() %>%
    filter(sample == samp) %>%
    select(v0, conc, sample) %>%
    do(fit_MM(data = .))
  MM.data <- bind_rows(MM.data, temp)
}

ggplot(MM.data, mapping = aes(x = conc, y = v0, color = sample)) + geom_point() + geom_line(aes(x = conc, y = predicted_v0))

