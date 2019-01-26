#### import libraries
library(tidyverse)
library(viridis)
library(lubridate)
calibration <- tibble("sensor_conc" = c(10, 20), 
                      "cal_slope" = c(3465, 2667), 
                      "cal_intercept" = c(3445, 5259))
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
    mutate("Time" = as.integer(as.numeric(hms(Time))))  #### hms gives a warning when parsing 00:00:00
  data <- index %>% 
    filter(., data_file == file) %>% 
    select(-data_file) %>% 
    inner_join(., data_gathered, by = "well") %>% 
    mutate("row" = str_sub(well, 1, 1), "column" = str_sub(well, 2)) %>% 
    mutate("condition" = str_c(date, sample, well, conc, GAP_conc, sensor_conc, sep = "-")) %>% 
    inner_join(., calibration, by = "sensor_conc") %>% 
    mutate("product_conc" = as.integer(fluorescence - cal_intercept)/cal_slope) %>% 
    filter(! is.na(fluorescence))
  return(data)
}

# combine all the files into one tibble
outfile <- "GTPase_assay/GAP_assay_data_parsed.txt"
index <- read_tsv("GTPase_assay/2018_data/data_index.txt", col_names = T)
files <- index %>% 
  pull(data_file) %>% unique()

### read in the data files, join them with the info from the index file and make them tidy
#data_points_to_discard <- read_tsv("GEF_assay/2018_data/data_to_discard.txt")
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
write_tsv(dataset, path = outfile)



run_linear <- function(data) {
  start_intercept <- data %>% pull(intercept) %>% unique()
  percent_fit <- data %>% pull(percent_fit) %>% unique()
  max_flo <- max(data$observed, na.rm = T)
  min_flo <- min(data$observed, na.rm = T)
  delta_flo <- max_flo - min_flo
  cutoff_fluo <- min_flo + (percent_fit * delta_flo)
  initial_data <- data %>% filter(observed < cutoff_fluo)
  data %>% ggplot(aes(x = Time, y = observed)) + geom_point() + geom_point(data = initial_data, color = "red")
  fit <- lm(I(observed - start_intercept) ~ Time - 1, data = initial_data)
  slope <- fit$coefficients
  data$predicted <- slope * data$Time + start_intercept  
  data$slope <- slope
  data$intercept <- start_intercept
  data$f0 <- start_intercept
  data$f_plateau <- max_flo
  data$vf0 <- - slope / (GAP_conc * 0.001)
  data$percent_fit <- percent_fit
  ## label the initial points used for fitting
  data <- data %>% 
    mutate("initial_linear" = ifelse(observed < cutoff_fluo, "initial_linear", "all")) %>% 
    mutate("predicted" = ifelse(predicted < max_flo, predicted, NA))
  ### calculate the initial rate
  data$vf0 <- - slope / (data$GAP_conc[1] * 0.001)  ### initial rate in fluorescence units
  return(data)
}
run_nls <- function(data) {
  data <- data %>% 
    filter(Time < cutoff_time)
  #max_flo <- max(data$observed, na.rm = T)
  #t_at_max <- data$Time[data$observed == max_flo]
  #data <- data %>% 
   # filter(Time < t_at_max)
  c0 <- data$conc[1]
  GAP_conc <- data$GAP_conc[1]
  start <- list(
    f0 = min(data$observed, na.rm = T),
    f_plat = max(data$observed, na.rm = T),
    k = 1e-2
    )
  lower <- c(
    min(0.7 * min(data$observed, na.rm = T), 1.2 * min(data$observed, na.rm = T)),
    0.8 * (max(data$observed, na.rm = T)),
    1e-5
    )
  upper <- c(
    max(0.7 * min(data$observed, na.rm = T), 1.2 * min(data$observed, na.rm = T)),
    1.3 * max(data$observed, na.rm = T),
    1
    )
  # f0 + (f_plat - f0) * (1 - exp(-k * Time))
  out <- nlsLM(observed ~ f0 + (f_plat - f0) * (1 - exp(-k * Time)),
               data = data,
               start = start,
               lower = lower,
               upper = upper,
               control = nls.lm.control(maxiter = 500))
  
  #### save optimal parameters
  f0 <- coef(out)[1]
  f_plat <- coef(out)[2]
  k <- coef(out)[3]
  print(str_c(data$condition[1], " ", k, " ", f0, " ", f_plat))
  k_background <- coef(out)[5]
  data$predicted <- f0 + (f_plat - f0) * (1 - exp(-k * data$Time))
  #### save optimal parameters in the data table
  data$k <- k
  data$f_plateau <- f_plat
  data$f0 <- f0 
  span <- f_plat - f0
  data$span <- span
  ## get the fit statistics
  # Chi Square in kaleidagraph is the total sum of the squared errors (sum((y-f(x))/sigma)^2)
  # R in kaleidagraph is the Pearson's R (root 1- chi^2/sum(sig*(y-mean(y))^2)
  sigma <- summary(out)$sigma
  resid <- summary(out)$residuals
  chisq <- out$m$deviance()
  data$chisq <- chisq
  data$pearsonr <- sqrt(1-chisq/sum(sigma*resid^2))
  data$span_pval <- log10(summary(out)$coefficients[,4][2])
  data$k_pval <- log10(summary(out)$coefficients[,4][4])
  ### calculate the initial rate
  data$vf0 = (span * k) / (GAP_conc * 0.001) ### initial rate in fluorescence units
  return(data)
}

plot_fits <- function(data, output = getwd()) {
  samples <- data %>% pull(sample) %>% unique()
  for (j in seq_along(samples)) {
    data_mutant <- data %>% filter(sample == samples[j])
    plots <- list()
    conditions <- data_mutant %>% 
      arrange(conc, condition) %>% pull(condition) %>% unique()
    for (i in seq_along(conditions)) {
      condition_data <- data_mutant %>% 
        filter(condition == conditions[i])
      fit <- condition_data %>% pull(fit) %>% unique()
      
      if (fit == "exp") {
        title <- condition_data %>% 
          select(condition, pearsonr, k, span, vf0) %>% 
          unique() %>%
          mutate("title" = str_c("exp fit\n", condition, "R =", round(pearsonr, 2), "k =", round(k, 4), "span =", round(span, 0), sep = " ")) %>% 
          pull(title)
      data_to_plot <- condition_data %>% 
        select(Time, observed, predicted) %>% 
        gather(`curve type`, value, -Time)
      plots[[i]] <- ggplot(data_to_plot[data_to_plot[["curve type"]] != "observed", ], mapping = aes(x = Time, y = value, color = `curve type`)) +
          geom_point(data = data_to_plot[data_to_plot[["curve type"]] == "observed", ], color = "black", alpha = 0.5) + 
          geom_line() + ylab("fluorescence") +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          ggtitle(title)
      } else if (fit == "lin") {
        
        title <- condition_data %>% 
          select(condition, slope, intercept) %>% 
          mutate("title" = str_c("linear fit\n", condition,  "slope =", round(slope, 1), sep = " ")) %>% 
          pull(title) %>% unique()
        slope <- condition_data %>% pull(slope) %>% unique()
        intercept <- condition_data %>% pull(intercept) %>% unique()
        data_to_plot <- condition_data %>% 
          select(Time, initial_linear, observed) %>% 
          mutate("data" = initial_linear)
        plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed, color = data)) +
          geom_point(alpha = 0.5) + 
          geom_abline(slope = slope, intercept = intercept) + ylab("fluorescence") +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          ggtitle(title)
        
      }
    }
    pdf(paste0(output, str_c(samples[j], 'fits.pdf', sep = '_')), width = 10)
    print(plots)
    dev.off()
  }
}

today <- gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
output <- str_c("GTPase_assay/", today, '_output', '/')
dir.create(output, showWarnings = FALSE)


exp_fits <- dataset %>% 
  filter(fit == "exp") %>% pull(condition) %>% unique()
if (length(exp_fits) > 0) {
  processed.data_exp <- dataset %>%
    filter(fit == "exp") %>% 
    group_by(sample, condition) %>%
    do(run_nls(.)) %>%  # fit curve
    ungroup() %>% 
    mutate("v0" = vf0/cal_slope)
}

lin_fits <- dataset %>% 
  filter(fit == "lin") %>% pull(sample) %>% unique()
if (length(lin_fits) > 0) {
  processed.data_lin <- dataset %>% 
    filter(fit == "lin") %>% 
    group_by(condition) %>%
    do(run_linear(.)) %>%  # linear fit
    ungroup() %>% 
    mutate("v0" = -1 * vf0/cal_slope)
}

if (length(lin_fits) > 0 & length(exp_fits) > 0) {
  processed.data <- full_join(processed.data_exp, processed.data_lin) 
} else if (length(exp_fits) > 0) {
  processed.data <- processed.data_exp
} else if (length(lin_fits) > 0) {
  processed.data <- processed.data_lin
}

processed.data %>% ggplot(aes(x = conc, y = v0, color = row, shape = as.character(sensor_conc))) + 
  geom_point()
#write_tsv(processed.data, file.path(output, "fit_data.txt"))
plot_fits(processed.data, output = output)




dataset %>% 
  filter(Time < 750 & sample == "PE1_WT") %>% 
  ggplot(aes(Time, observed, color = conc)) + 
  geom_point() + facet_wrap(~sensor_conc) +
  scale_color_viridis()
