#### import libraries
library(tidyverse)
library(lubridate)
library(minpack.lm)

run_nls <- function(data, deadtime, debug = FALSE) {
  
  c0 <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  
  #### debug mode: if debug is TRUE, print out condition being fit 
  if (debug) { print(unique(data$condition)) }

  start <- list(f_plateau = min(data$observed, na.rm = T),
                span1 = max(data$observed, na.rm = T) - min(data$observed, na.rm = T),
                span2 = 0.1*(max(data$observed, na.rm = T) - min(data$observed, na.rm = T)),
                k = 5e-3,
                k_background = 0)
  lower <- c(min(data$observed, na.rm = T),
             0.0,
             0.0,
             0.0,
             0.0)
  upper <- c(max(data$observed, na.rm = T),
             max(data$observed, na.rm = T),
             max(data$observed, na.rm = T),
             1.0,
             1.0)
  
  out <- nlsLM(observed ~ span1 * exp(-k * (Time+deadtime)) - span2 * (1 - exp(-k_background * (Time+deadtime))) + f_plateau,
               data = data,
               start = start,
               lower = lower,
               upper = upper,
               control = nls.lm.control(maxiter = 300))
  
  #### save optimal parameters
  f_plateau <- coef(out)[1]
  span1 <- coef(out)[2]
  span2 <- coef(out)[3]
  k <- coef(out)[4]
  k_background <- coef(out)[5]

  data$predicted <- span1 * exp(-k * (data$Time+deadtime)) - span2 * (1 - exp(-k_background * (data$Time+deadtime))) + f_plateau
  data$exchange <- span1 * exp(-k * (data$Time+deadtime))  + f_plateau
  data$background <-  - span2 *(1 - exp(-k_background * (data$Time+deadtime)))

  #### save optimal parameters in the data table
  data$k <- k
  data$f_plateau <- f_plateau
  data$k_background <- k_background
  data$span1 <- span1
  data$span2 <- span2

  ### calculate the initial rate
  data$vf0 = (span1 * k * exp(k * 0)) / (GEF_conc*0.001)   ### initial rate in fluorescence units

  return(data)
}

fit_conversion_factor <- function(data) {
  c_plateau_values <- list()
  f_plateau_values <- list()
  conditions = unique(data$condition)
  for (i in seq_along(conditions)) {
    data_temp <- data %>% filter(condition == conditions[i])
    c_plateau_values[[i]] = 0.995*data_temp$conc[1]  # account 1/200 molecules being GDP instead of mGTP at steady-state
    f_plateau_values[[i]] = data_temp$f_plateau[1]
  }
  c_plateau_values <- unlist(c_plateau_values)
  f_plateau_values <- unlist(f_plateau_values)

  out <- lm(c_plateau_values ~ f_plateau_values)
  data$conversion_ratio = out$coef[2] # need to divide by the factor that accounts for trp signal decrease for mant bound
  data$v0 = data$vf0 * data$conversion_ratio

  return(data)
}

fit_MM <- function(data) {
  data_to_fit <- data %>% filter(Time == 0)
  start <- list(Vmax = 5, Km = 2)
  out <- nlsLM(v0 ~ (conc * Vmax) / (Km + conc),
               data = data_to_fit,
               start = start,
               control = nls.lm.control(maxiter = 200))
  Vmax <- coef(out)[1]
  Km <- coef(out)[2]
  data$kcat <- Vmax
  data$Km <- Km
  data$predicted_v0 <- (data$conc * Vmax) / (Km + data$conc)
  return(data)
}

plot_raw_data <- function(data, by_sample = FALSE, output = getwd()) {
  
  today = gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
  
  if (by_sample) {
    samples = unique(data$sample)
    for (j in seq_along(samples)) {
        data_mutant <- data %>% filter(sample == samples[j])
        plots <- list()
        conditions = unique(data_mutant$condition)
        for (i in seq_along(conditions)) {
          data_to_plot <- data_mutant %>% filter(condition == conditions[i]) 
          plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
            geom_point(color = "black") +
            ggtitle(conditions[i])
        }
        pdf(paste0(output, paste(samples[j], 'raw_data.pdf', sep = '_')))
        print(plots)
        dev.off()
    }
  if (!by_sample) {
    plots <- list()
    conditions = unique(data$condition)
    for (i in seq_along(conditions)) {
      data_to_plot <- data %>% filter(condition == conditions[i]) 
      plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
        geom_point(color = "black") +
        ggtitle(conditions[i])
    }
    pdf(paste0(output, 'raw_data.pdf'))
    print(plots)
    dev.off()
    }
  }
}

plot_fits <- function(data, by_sample = FALSE, output = getwd()) {
  
  if (by_sample) {
    samples = unique(data$sample)
    for (j in seq_along(samples)) {
        data_mutant <- data %>% filter(sample == samples[j])
        plots <- list()
        conditions = unique(data_mutant$condition)
        for (i in seq_along(conditions)) {
          data_to_plot <- data_mutant %>% filter(condition == conditions[i]) 
          plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
            geom_point(color = "black") +
            geom_line(data_to_plot, mapping = aes(x = Time, y = predicted), color = "red") +
            geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
            ggtitle(conditions[i])
        }
        pdf(paste0(output, paste(samples[j], 'fits.pdf', sep = '_')))
        print(plots)
        dev.off()
    }
  if (!by_sample) {
    plots <- list()
    conditions = unique(data$condition)
    for (i in seq_along(conditions)) {
      data_to_plot <- data %>% filter(condition == conditions[i]) 
      plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
        geom_point(color = "black") +
        geom_line(data_to_plot, mapping = aes(x = Time, y = predicted), color = "red") +
        geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
        ggtitle(conditions[i])
    }
    pdf(paste0(output, 'fits.pdf'))
    print(plots)
    dev.off()
    }
  }
}

plot_fits_show_bkgrnd <- function(data, by_sample = FALSE, output = getwd()) {
  if (by_sample) {
    samples = unique(data$sample)
    for (j in seq_along(samples)) {
        data_mutant <- data %>% filter(sample == samples[j])
        plots <- list()
        conditions = unique(data_mutant$condition)
        for (i in seq_along(conditions)) {
          data_to_plot <- data_mutant %>% filter(condition == conditions[i]) 
          plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
            geom_point(color = "black") +
            geom_line(data_to_plot, mapping = aes(x = Time, y = predicted), color = "red") +
            geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
            geom_line(data_to_plot, mapping = aes(x = Time, y = background), color = "blue") +
            ggtitle(conditions[i])
        }
        pdf(paste0(output, paste(samples[j], 'fits_show_bkgrnd.pdf', sep = '_')))
        print(plots)
        dev.off()
    }
  if (!by_sample) {
    plots <- list()
    conditions = unique(data$condition)
    for (i in seq_along(conditions)) {
      data_to_plot <- data %>% filter(condition == conditions[i]) 
      plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
        geom_point(color = "black") +
        geom_line(data_to_plot, mapping = aes(x = Time, y = predicted), color = "red") +
        geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
        geom_line(data_to_plot, mapping = aes(x = Time, y = background), color = "blue") +
        ggtitle(conditions[i])
    }
    pdf(paste0(output, 'fits_show_bkgrnd.pdf'))
    print(plots)
    dev.off()
    }
  }
}

plot_MM <- function(data, output = getwd()) {
  plots <- list()
  samples = unique(data$sample)
  for (i in seq_along(samples)) {
    data_to_plot <- data %>% filter(sample == samples[i]) %>% filter(Time == 0)
    plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = conc, y = v0, color = as.character(date))) + 
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0, color = sample)) +
      ggtitle(paste(samples[i], paste("Km:", data_to_plot$Km[1], sep = " "), paste("kcat:", data_to_plot$kcat[1], sep = " "), sep = "\n"))
  }
  pdf(paste0(output, 'MM.pdf'))
  print(plots)
  dev.off()
}

plot_parameters <- function(data, output = getwd()) {
  samples = unique(data$sample)
  for (i in seq_along(samples)) {
    data_to_plot <- data %>% filter(sample == samples[i]) %>% filter(Time == 0)
    pdf(paste0(output, paste(samples[i], 'params.pdf', sep = '_')))
    print(ggplot(data_to_plot, aes(x = conc, y = 1000*k, color = as.character(date))) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = f_plateau, color = as.character(date))) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = k_background, color = as.character(date))) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = span1, color = as.character(date))) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = span2, color = sample)) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = vf0, color = as.character(date))) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = v0, color = as.character(date))) + geom_point())
    print(ggplot(data_to_plot, aes(x = conc, y = v0/0.4, color = as.character(date))) + geom_point())
    # print(ggplot(data_to_plot, aes(x = conc, y = v0, color = as.character(date))) + geom_point() + geom_text(aes(label=conc),hjust=0, vjust=0) + coord_cartesian(xlim = c(0,15), ylim = c(0, 6)))
    dev.off()
  }
}

#### set output directory
today = gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
output <- paste0(today, '_output', '/')
dir.create(file.path(getwd(), output), showWarnings = FALSE)

#### read in data
inputfilename <- "GEF_assay/biotek_test_exp/data/good_data_parsed.txt"
biotek.data <- read_delim(inputfilename, delim = "\t", col_names = T)

#### cut off curves with odd bumps
cutoff_curves <- c('20180213-PE9_D79S-A5-3.75',
                '20180213-PE9_D79S-A6-4',
                '20180215-PE1_WT-A4-3',
                '20180215-PE1_WT-A5-3.75',
                '20180206-PE1_WT-E5-2.25',
                '20180215-PE9_D79S-A8-2')
biotek.data <- subset(biotek.data, ! condition %in% cutoff_curves | Time < cutoff_time)

#### Fit data assuming photobleaching decay is exponential, with observations in fluorescence units
processed.data <- biotek.data %>%
  group_by(condition) %>%
  mutate(observed = fluorescence) %>%
  do(run_nls(., deadtime = 0, debug = F)) %>%  # fit curve
  ungroup() %>%
  group_by(date, row) %>%
  do(fit_conversion_factor(.)) %>%  # fit conversion ration from f_plateau
  ungroup() %>%
  group_by(sample) %>%
  do(fit_MM(.))  # fit michaelis-menten

# plot_raw_data(processed.data, by_sample = T, output = output)
# plot_fits(processed.data, by_sample = T, output = output)
# plot_fits_show_bkgrnd(processed.data, by_sample = T, output = output)
# plot_MM(processed.data, output = output)
# plot_parameters(processed.data, output = output)
