#### import libraries
library(tidyverse)
library(lubridate)
library(minpack.lm)

remove(list = ls())

## READ and FORMAT the raw ASCII data from Synergy H1
####### a function to read in all the biotek files, gather the data and reformat the time columns
read_and_gather <- function(file) {
  print(file)
  raw_in <-  read_lines(file)
  first_data_row <- grep(raw_in, pattern = "Time\tT")
  last_data_row <- grep(raw_in, pattern = "Results")
  if (length(last_data_row) > 0) {
    n_rows_to_read <- last_data_row - first_data_row - 2
    data_in <-  read_tsv(file, col_names = T, skip = (first_data_row - 1),
                           n_max = n_rows_to_read, locale = locale(encoding = 'windows-1250'))
  } else {
    data_in <-  read_tsv(file, col_names = T, skip = (first_data_row - 1),
                           locale = locale(encoding = 'UTF-8'))
  }
  data_gathered <- data_in %>% 
    select(., -2) %>% 
    gather(., key = well, value = fluorescence, -Time) %>% 
    mutate("Time" = as.numeric(hms(Time)))  #### hms gives a warning when parsing 00:00:00
  dataset <- index %>% 
    filter(., data_file == file) %>% 
    inner_join(., data_gathered, by = "well") %>% 
    #select(., -data_file) %>% 
    mutate("condition" = str_c(date, sample, well, conc, sep = "-"),
           "row" = str_sub(well, 1, 1), "column" = str_sub(well, 2)) %>% 
    mutate("cutoff_time" = ifelse(is.na(cutoff_time), max(Time, na.rm = T), cutoff_time)) %>% 
    filter(Time < cutoff_time)
  return(dataset)
}

# combine all the files into one tibble
#(files <- dir("GEF_assay/2018_data", pattern = "GEF_FRET_assay", full.names = T))
outfile <- "GEF_assay/good_data_parsed.txt"
### load the index file (has conditions per well)
(index <- read_tsv("GEF_assay/2018_data/data_index.txt", col_names = T))
#exemplary_data <- read_tsv("GEF_assay/2018_data/exemplary_data.txt", col_names = F) %>% 
 # pull(X1)
files <- index %>% 
  pull(data_file) %>% unique()

### read in the data files, join them with the info from the index file and make them tidy
#data_points_to_discard <- read_tsv("GEF_assay/2018_data/data_to_discard.txt")
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
write_tsv(dataset, path = outfile)
#dataset <- read_tsv(outfile)
### plot fluorescence versus concentration calibration
# calibration_data <- dataset %>% 
#   filter(date == "20181003" & k_est == 0)
# calibration_conc <- calibration_data %>% pull(conc) %>% unique()
# back_k_calibration_table <- data.frame()
# plots <- list()
# for (i in seq_along(calibration_conc)) {
#   cc <- calibration_conc[i]
#   data <- calibration_data %>% filter(conc == cc)
#   plots[[i]] <- data %>% ggplot(aes(x = Time, y = fluorescence)) + geom_point()
#   max_flo <- max(data$fluorescence, na.rm = T)
#   start <- list(f_plateau = 0.3 * max(data$fluorescence, na.rm = T),
#                 span = max(data$fluorescence, na.rm = T) - 0.2 * max(data$fluorescence, na.rm = T),
#                 k = 1e-4)
#   lower <- c(f_plateau = 0.2 * max(data$fluorescence, na.rm = T),
#              span = max(data$fluorescence, na.rm = T) - 0.7 * max(data$fluorescence, na.rm = T),
#              k = 1e-6)
#   upper <- c(f_plateau = 0.7 * max(data$fluorescence, na.rm = T),
#              span = max(data$fluorescence, na.rm = T) - 0.2 * max(data$fluorescence, na.rm = T),
#              k = 1e-3)
#   data %>% ggplot(aes(x = Time, y = fluorescence)) + geom_point()
#   out <- nlsLM(fluorescence ~ span * exp(-k * (Time)) + f_plateau,
#                    data = data, control = nls.lm.control(maxiter = 500), 
#                    start = start, lower = lower, upper = upper) %>% 
#   print()
#   back_k_calibration_table <- rbind(back_k_calibration_table, 
#             data.frame("conc" = cc, "f_plateau" = coef(out)[1], "span" = coef(out)[2], "k" = coef(out)[3], 
#                        "max" = max_flo, "span/max" = coef(out)[2]/max_flo
#                        ))
# }
# back_k_calibration_table <- as_tibble(back_k_calibration_table) 
# back_k_calibration_table %>% 
#   ggplot(aes(x = conc, y = f_plateau)) + geom_point()
# range(back_k_calibration_table$k)
dataset <- dataset %>% 
  filter(k_est != 0)

run_nls <- function(data, debug = FALSE) {
  
  c0 <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  k_est <- data$k_est[1]
  span1_est <- data$span1[1]
  span2_est <- 283.5 * c0
  f_plateau_est <- 1849.6  * c0
  f_plateau_est <- data$f_plateau[1]
  #### debug mode: if debug is TRUE, print out condition being fit 
  if (debug) { print(unique(data$condition)) }
  start <- list(
                f_plateau = f_plateau_est,
                span1 = span1_est,
                span2 = span2_est,
                k = k_est,
                k_background = 1e-4)
  lower <- c(
                f_plateau_est - 0.2 * f_plateau_est,
                span1_est - 0.1*span1_est,
               0,
              1e-5,
               1e-5)
  upper <- c(
              f_plateau_est + 0.2 * f_plateau_est,
             span1_est + 0.1*span1_est,
              max(data$observed, na.rm = T),
              0.1,
               3e-4)

  out <- nlsLM(observed ~ span1 * exp(-k * Time) + span2 * (exp(-k_background * Time)) + f_plateau,
               data = data,
               start = start,
               lower = lower,
               upper = upper,
               control = nls.lm.control(maxiter = 500))
  
  #### save optimal parameters
  f_plateau <- coef(out)[1]
  span1 <- coef(out)[2]
  span2 <- coef(out)[3]
  k <- coef(out)[4]
  print(str_c(data$condition[1], " ", k, " ", f_plateau, " ", span1))
  k_background <- coef(out)[5]

  data$predicted <- span1 * exp(-k * data$Time) + span2 * (exp(-k_background * data$Time)) + f_plateau
  
  data$exchange <- span1 * exp(-k * data$Time)  + f_plateau + span2
  data$background <-  span2 * (exp(-k_background * data$Time)) + f_plateau #+ span1

  #### save optimal parameters in the data table
  data$max_flo <- max(data$observed, na.rm = T)
  data$k <- k
  data$f_plateau <- f_plateau
  data$k_background <- k_background
  data$span1 <- span1
  data$span2 <- span2
  f0 <- span1 + span2 + f_plateau
  data$f_mid <- f0 - span1
  data$f0 <- f0

  ### calculate the initial rate
  data$vf0 <- (span1 * k * exp(k * 0)) / (GEF_conc*0.001)   ### initial rate in fluorescence units
  return(data)
}

fit_conversion_factor <- function(data) {
  sample <- data$sample[1]
  date <- data$date[1]
  c_plateau_values <- list()
  f_plateau_values <- list()
  conditions = unique(data$condition)
  for (i in seq_along(conditions)) {
    data_temp <- data %>% filter(condition == conditions[i])
    c_plateau_values[[i]] = 0.995*data_temp$conc[1]  # account 1/200 molecules being GDP instead of mGTP at steady-state
    f_plateau_values[[i]] = data_temp$f_mid[1]
  }
  c_plateau_values <- unlist(c_plateau_values)
  f_plateau_values <- unlist(f_plateau_values)
  tib <- tibble("conc" = c_plateau_values, "fluor" = f_plateau_values)
  out <- lm(c_plateau_values ~ f_plateau_values)
  tib %>% ggplot(aes(fluor, conc)) + geom_point() + 
    geom_abline(slope = out$coeff[2], intercept = out$coeff[1]) +
    ggtitle(str_c(sample, " ", date))
  data$conversion_ratio = out$coef[2] # need to divide by the factor that accounts for trp signal decrease for mant bound
  data$v0 = data$vf0 * data$conversion_ratio
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
          plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = fluorescence)) +
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
            # geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
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
        # geom_line(data_to_plot, mapping = aes(x = Time, y = exchange), color = "green") +
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
    samples <- data %>% pull(sample) %>% unique()
    for (j in seq_along(samples)) {
        data_mutant <- data %>% filter(sample == samples[j])
        plots <- list()
        conditions <- data_mutant %>% 
          arrange(conc, condition) %>% pull(condition) %>% unique()
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
    conditions <- data %>% 
      arrange(conc, condition) %>% pull(condition) %>% unique()
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
  sample <- data$sample[1]
  print(sample)
  plots <- list()
  samples = unique(data$sample)
  for (i in seq_along(samples)) {
    data_to_plot <- data %>% filter(sample == samples[i]) %>% filter(Time == 0)
    plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = conc, y = v0, color = as.character(date))) + 
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0, color = sample)) +
      ggtitle(paste(samples[i], paste("Km:", round(data_to_plot$Km[1],2), sep = " "), paste("kcat:", round(data_to_plot$kcat[1],2), sep = " "), sep = "\n"))
  }
  pdf(paste0(output, 'MM.pdf'), width = 10)
  print(plots)
  dev.off()
}

plot_MM_bins <- function(data, output = getwd()) {
  
  plots <- list()
  samples = unique(data$sample)
  merged_data <- data.frame()
  for (i in seq_along(samples)) {
    
    data_to_plot <- data %>%
      filter(sample == samples[i]) %>%
      filter(Time == 0) %>% 
      mutate("floor_conc" = ifelse(conc > 1, floor(conc), conc))
    
   data_to_plot <- data_to_plot %>%
      group_by(floor_conc) %>%
      summarise(mean_v0 = mean(v0)) %>%
      inner_join(data_to_plot, by = "floor_conc")

    data_to_plot <- data_to_plot %>%
      group_by(floor_conc) %>%
      summarise(sd_v0 = sd(v0)) %>%
      inner_join(data_to_plot, by = "floor_conc")
    merged_data <- rbind(merged_data, data_to_plot)

    plots[[i]] <- ggplot(data_to_plot, aes(x = floor_conc, y = mean_v0)) +
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0)) +
      geom_errorbar(aes(ymin = mean_v0 - sd_v0, ymax = mean_v0 + sd_v0)) +
      ggtitle(paste(samples[i],
                paste("Km:", round(data_to_plot$Km[1],2), sep = " "),
                paste("kcat:", round(data_to_plot$kcat[1],2), sep = " "),
              sep = "\n"))
  }
  plots[[i + 1]] <- ggplot(merged_data, aes(x = floor_conc, y = mean_v0, color = sample)) +
    geom_point() +
    geom_line(aes(x = conc, y = predicted_v0, color = sample)) +
    geom_errorbar(aes(ymin = mean_v0 - sd_v0, ymax = mean_v0 + sd_v0, color = sample))
  pdf(paste0(output, 'MM_bins.pdf'), width = 10)
  print(plots)
  dev.off()
  
  for (i in seq_along(samples)) {
    
    data_to_plot <- data %>%
      filter(sample == samples[i]) %>%
      filter(Time == 0)
    
   data_to_plot <- data_to_plot %>%
      group_by(conc) %>%
      summarise(mean_v0 = mean(v0)) %>%
      inner_join(data_to_plot, by = "conc")

    data_to_plot <- data_to_plot %>%
      group_by(conc) %>%
      summarise(sd_v0 = sd(v0)) %>%
      inner_join(data_to_plot, by = "conc")

    plots[[i]] <- ggplot(data_to_plot, aes(x = conc, y = mean_v0)) +
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0, color = sample)) +
      ggtitle(paste(samples[i],
                paste("Km:", round(data_to_plot$Km[1],2), sep = " "),
                paste("kcat:", round(data_to_plot$kcat[1],2), sep = " "),
              sep = "\n"))
  }
  pdf(paste0(output, 'MM_bins_no_error_bars.pdf'))
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
today <- gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
output <- str_c(today, '_output', '/')
dir.create(output, showWarnings = FALSE)

### first plot raw data for everything
#plot_raw_data(dataset, by_sample = T, output = output)

# dataset <- dataset %>% 
#   filter( ! condition %in% data_points_to_discard$condition)
#   #filter(condition %in% exemplary_data)
# 

#### Fit data assuming photobleaching decay is exponential, with observations in fluorescence units
processed.data <- dataset %>%
  #filter(sample == "PE15_T34L") %>% 
  group_by(sample, condition) %>%
  mutate(observed = fluorescence) %>%
  do(run_nls(., debug = T)) %>%  # fit curve
  ungroup() %>%
  group_by(sample, date) %>% 
  #group_by(sample, row, date) %>%
  do(fit_conversion_factor(.)) %>%  # fit conversion ratio from f_plateau
  ungroup()


### make a function where I can quickly normalize and plot raw data for a couple of selected measurements
## use it when something can't be fit, to compare to the curves for similar concentrations...
# plot_selected_fits <- function(sam, min_conc, max_conc) {
#   to_plot <- processed.data %>% 
#     filter(sample == sam & conc > min_conc & conc < max_conc)
#   to_plot %>% ggplot(aes(x = Time, y = observed, col = as.character(condition))) + 
#     geom_point() +
#     geom_line(aes(x = Time, y = exchange))
# }
# plot_selected_fits("PE1_WT", 11, 13)
# plot_norm_selected_fits <- function(sam, min_conc, max_conc) {
#   to_plot <- processed.data %>%
#     filter(sample %in% sam & conc > min_conc & conc < max_conc & Time < 1500) %>%
#     group_by(condition) %>%
#     mutate("norm_fluorescence" = 1 - ((quantile(fluorescence, prob = 0.999, na.rm = T) - fluorescence) /
#                                         (quantile(fluorescence, prob = 0.999, na.rm = T) - quantile(fluorescence, prob = 0.001, na.rm = T))) ) %>%
#     ungroup()
#   to_plot %>% ggplot(aes(x = Time, y = norm_fluorescence, col = as.character(condition))) +
#     geom_point()
# }
# plot_norm_selected_fits(c("PE2_T34A"), 2.5, 3.5)

#plot_fits(processed.data, by_sample = T, output = output)
#plot_fits_show_bkgrnd(processed.data, by_sample = T, output = output)

fit.parameters <- processed.data %>% 
  #select(sample, date, well, conc, k, f_plateau, k_background, span1, span2, vf0, conversion_ratio, v0) %>% 
  select(sample, condition, date, well, conc, k, f_plateau, f_mid, f0,  k_background, span1, span2, max_flo, vf0, v0) %>% 
  unique() %>% 
  arrange(sample, conc)

#write_tsv(fit.parameters, file.path(output, "fit_parameters.txt"))
fit.parameters %>%
  
  #filter(sample == "PE1_WT") %>%
  ggplot(., aes(x = conc, y = v0, color = as.character(date))) + geom_point() + facet_wrap(~ sample)

# fit <- lm(fit.parameters$span1 ~ fit.parameters$conc)
# fit.parameters %>%
#    ggplot(., aes(x = conc, y = span1, color = as.character(date), shape = sample)) +
#    geom_point() + geom_abline(aes(slope = fit$coefficients[2], intercept = fit$coefficients[1]))
# #
# #
#  fit <- lm(fit.parameters$f_plateau ~ fit.parameters$conc)
#  fit.parameters %>%
#    #filter(sample == "PE1_WT") %>%
#    ggplot(., aes(x = conc, y = f_plateau, color = as.character(date), shape = sample)) +
#    geom_point() + geom_abline(aes(slope = fit$coefficients[2], intercept = fit$coefficients[1]))

# fit <- lm(fit.parameters$span2 ~ fit.parameters$k_background)
# fit.parameters %>%
#   #filter(sample == "PE1_WT") %>%
#   ggplot(., aes(x = k_background, y = span2, color = as.character(date), shape = sample)) + 
#   geom_point() + geom_abline(aes(slope = fit$coefficients[2], intercept = fit$coefficients[1]))

# fit.parameters %>% 
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = conc, y = k_background, color = sample)) + geom_point()
# fit.parameters %>%   
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = conc, y = k, color = as.character(date))) + geom_point() + facet_wrap(~ sample)
# fit.parameters %>%   
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = conc, y = span1, color = as.character(date))) + geom_point() + facet_wrap(~ sample)
# 
# fit.parameters %>%   
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = conc, y = vf0, color = as.character(date))) + geom_point() + facet_wrap(~ sample)
# 
# fit.parameters %>%   
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = conc, y = f_plateau, color = as.character(date))) + geom_point() + facet_wrap(~ sample)
# 
# fit.parameters %>%   
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = conc, y = v0, color = as.character(date))) + geom_point() + facet_wrap(~ sample)
# 
# fit.parameters %>%   
#   filter(sample == "PE1_WT") %>% 
#   ggplot(., aes(x = max_flo, y = span1, color = as.character(date))) + geom_point() + facet_wrap(~ sample)
# 
# 
# 
fit_MM <- function(data) {
  sample <- data$sample[1]
  print(sample)
  data_to_fit <- data %>% filter(Time == 0)
  #start <- list(Vmax = 5, Km = 2)
  start <- list(Vmax = max(data_to_fit$v0, na.rm = T), Km = 2)
  lower <- c(0.5 * max(data_to_fit$v0, na.rm = T), 0.5)
  upper <- c(5 * max(data_to_fit$v0), 20)
  data_to_fit %>% ggplot(., aes(conc, v0)) + geom_point()
  out <- nlsLM(v0 ~ (conc * Vmax) / (Km + conc),
               data = data_to_fit,
               start = start,
               lower = lower,
               upper = upper,
               control = nls.lm.control(maxiter = 1000))
  Vmax <- coef(out)[1]
  Km <- coef(out)[2]
  data$kcat <- Vmax
  data$Km <- Km
  data$kcat_sd <- summary(out)$parameters[1,2]
  data$Km_sd <- summary(out)$parameters[2,2]
  data$predicted_v0 <- (data$conc * Vmax) / (Km + data$conc)
  return(data)
}

MM.data <- processed.data %>%
  #filter(sample != "PE12_H141R" & sample != "PE14_Y157A" & sample != "PE6_H141V") %>%
  group_by(sample) %>%
  do(fit_MM(.))  # fit michaelis-menten


#plot_MM(MM.data, output = output)

MM.data %>% ggplot(., aes(conc, v0, color = sample)) + geom_point() +
  geom_line(aes(x = conc, y = predicted_v0, color = sample))
# MM.data %>%
#   filter(sample %in% c("PE1_WT", "PE2_T34A", "PE3_T34E", "PE15_T34L", "PE11_T34G")) %>%
#   ggplot(., aes(conc, v0, color = sample)) + geom_point() +
#   geom_line(aes(x = conc, y = predicted_v0, color = sample))

#plot_parameters(processed.data, output = output)
#plot_MM_bins(MM.data, output = output)

MM.data %>% 
  select(sample, kcat, Km, kcat_sd, Km_sd) %>%
  unique() %>% 
  mutate() %>% 
  ggplot(aes(sample, kcat)) + geom_bar(stat = "identity") + geom_errorbar(aes(ymin = kcat - kcat_sd, ymax = kcat + kcat_sd))

MM.data %>% 
  select(sample, kcat, Km, kcat_sd, Km_sd) %>%
  unique() %>% 
  mutate() %>% 
  ggplot(aes(sample, Km)) + geom_bar(stat = "identity") + geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd))

# 
# #### analyse the control for interval time
# # check_time_interval <- dataset %>%
# #   filter(date == "20180906" & (row == "A" | row == "C") ) %>%
# #   filter(Time < 8000) %>%
# #   group_by(condition) %>%
# #   mutate("norm_fluorescence" = 1 - ((quantile(fluorescence, prob = 0.999, na.rm = T) - fluorescence) /
# #                                       (quantile(fluorescence, prob = 0.999, na.rm = T) - quantile(fluorescence, prob = 0.001, na.rm = T))) ) %>%
# #   ungroup()
# # 
# # check_time_interval %>%
# #   ggplot(., aes(x = Time, y = norm_fluorescence, color = row)) +
# #   geom_point() + facet_grid(sample ~ conc)
# # ggsave(str_c(output, "check_time_interval_norm_fluor.pdf"), width = 12)
# # check_time_interval %>%
# #   ggplot(., aes(x = Time, y = fluorescence, color = row)) +
# #   geom_point() + facet_grid(sample ~ conc)
# # ggsave(str_c(output, "check_time_interval.pdf"), width = 12)