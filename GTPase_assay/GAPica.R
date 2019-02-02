#### import libraries
library(tidyverse)
library(viridis)
library(lubridate)
calibration <- tibble("sensor_conc" = c(10, 20), 
                      "cal_slope" = c(0.0002841, 0.0003713), 
                      "cal_intercept" = c(-0.0014423, 0.1462813))
deadtime <- 30
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
    mutate("condition" = str_c(date, sample, well, conc, GAP_conc, sensor_conc, sep = "-")) 
  buffer_blanks <- data %>% filter(conc == 0) %>% 
    select(sensor_conc, Time, fluorescence)
  data <- data %>% filter(conc != 0) %>% 
    inner_join(., buffer_blanks, by = c("sensor_conc", "Time")) %>% 
    mutate("fluorescence" = fluorescence.x - fluorescence.y) %>% 
    rename("raw_fluo" = fluorescence.x, "background_fluo" = fluorescence.y) %>% 
    inner_join(., calibration, by = "sensor_conc") %>% 
    mutate("product_conc" = as.integer(fluorescence)*cal_slope + cal_intercept) %>% 
    filter(! is.na(fluorescence)) %>% 
    mutate("Time" = Time + deadtime)
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

dataset %>% 
  filter( (conc > 5 & conc < 10) | conc == 0) %>% 
  ggplot(aes(Time, fluorescence, color = as.character(conc))) + 
  geom_point()

write_tsv(dataset, path = outfile)


run_linear <- function(data) {
  condition <- data %>% pull(condition) %>% unique()
  print(condition)
  GAP_conc <- data$GAP_conc[1]
  start_intercept <- data %>% pull(intercept) %>% unique()
  percent_fit <- data %>% pull(percent_fit) %>% unique()
  max_flo <- max(data$product_conc, na.rm = T)
  min_flo <- min(data$product_conc, na.rm = T)
  delta_flo <- max_flo - min_flo
  cutoff_fluo <- min_flo + (percent_fit * delta_flo)
  initial_data <- data %>% filter(product_conc < cutoff_fluo)
  data %>% ggplot(aes(x = Time, y = product_conc)) + geom_point() + geom_point(data = initial_data, color = "red")
  fit <- lm(I(product_conc - start_intercept) ~ Time - 1, data = initial_data)
  slope <- fit$coefficients
  data$predicted <- slope * data$Time + start_intercept  
  data$slope <- slope
  data$intercept <- start_intercept
  data$f0 <- start_intercept
  data$f_plateau <- max_flo
  data$v0 <- slope / (GAP_conc * 0.001)
  data$percent_fit <- percent_fit
  ## label the initial points used for fitting
  data <- data %>% 
    mutate("initial_linear" = ifelse(product_conc < cutoff_fluo, "initial_linear", "all")) %>% 
    mutate("predicted" = ifelse(predicted < max_flo, predicted, NA))
  ### calculate the initial rate
  return(data)
}
run_nls <- function(data) {
  condition <- data %>% pull(condition) %>% unique()
  print(condition)
  cutoff_time <- data %>% pull(cutoff_time) %>% unique()
  if(! is.na(cutoff_time)) {
    data <- data %>% 
      filter(Time < cutoff_time)
  }
  c0 <- data$conc[1]
  GAP_conc <- data$GAP_conc[1]
  #### starting parameter estimates
  f0_start <- data %>% pull(intercept) %>% unique()
  f_plat_start <- data %>% pull(f_plateau_est) %>% unique()
  k_start <- data %>% pull(k_est) %>% unique()
  start <- list(
    f0 = f0_start,
    f_plat = f_plat_start,
    k = k_start
    )
  lower <- c(
    0.9 * f0_start,
    0.9 * f_plat_start,
    0.9 * k_start
    )
  upper <- c(
    1.1 * f0_start,
    1.1 * f_plat_start,
    1.1 * k_start
    )
  # f0 + (f_plat - f0) * (1 - exp(-k * Time))
  out <- nlsLM(product_conc ~ f0 + (f_plat - f0) * (1 - exp(-k * Time)),
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
  data$v0 = (span * k) / (GAP_conc * 0.001) ### initial rate in fluorescence units
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
          select(condition, pearsonr, k, span, v0) %>% 
          unique() %>%
          mutate("title" = str_c("exp fit\n", condition, "k =", round(k, 4), "span =", round(span, 3), sep = " ")) %>% 
          pull(title)
      data_to_plot <- condition_data %>% 
        select(Time, product_conc, predicted) %>% 
        gather(`curve type`, value, -Time)
      plots[[i]] <- ggplot(data_to_plot[data_to_plot[["curve type"]] != "product_conc", ], mapping = aes(x = Time, y = value, color = `curve type`)) +
          geom_point(data = data_to_plot[data_to_plot[["curve type"]] == "product_conc", ], color = "black", alpha = 0.5) + 
          geom_line() + ylab("product (free Pi) concentration") +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          ggtitle(title)
      } else if (fit == "lin") {
        
        title <- condition_data %>% 
          select(condition, slope, intercept) %>% 
          mutate("title" = str_c("linear fit\n", condition,  "slope =", round(slope, 4), sep = " ")) %>% 
          pull(title) %>% unique()
        slope <- condition_data %>% pull(slope) %>% unique()
        intercept <- condition_data %>% pull(intercept) %>% unique()
        data_to_plot <- condition_data %>% 
          select(Time, initial_linear, product_conc) %>% 
          mutate("data" = initial_linear)
        plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = product_conc, color = data)) +
          geom_point(alpha = 0.5) + 
          geom_abline(slope = slope, intercept = intercept) + ylab("product (free Pi) concentration") +
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
output <- str_c("GTPase_assay/", today, 'GAP_output', '/')
dir.create(output, showWarnings = FALSE)


exp_fits <- dataset %>% 
  filter(fit == "exp") %>% pull(condition) %>% unique()
if (length(exp_fits) > 0) {
  processed.data_exp <- dataset %>%
    filter(fit == "exp") %>% 
    group_by(sample, condition) %>%
    do(run_nls(.)) %>%  # fit curve
    ungroup() 
}

lin_fits <- dataset %>% 
  filter(fit == "lin") %>% pull(sample) %>% unique()
if (length(lin_fits) > 0) {
  processed.data_lin <- dataset %>% 
    filter(fit == "lin") %>% 
    group_by(condition) %>%
    do(run_linear(.)) %>%  # linear fit
    ungroup()
}

if (length(lin_fits) > 0 & length(exp_fits) > 0) {
  processed.data <- full_join(processed.data_exp, processed.data_lin) 
} else if (length(exp_fits) > 0) {
  processed.data <- processed.data_exp
} else if (length(lin_fits) > 0) {
  processed.data <- processed.data_lin
}

processed.data %>% 
  select(condition, conc, v0, row, sensor_conc) %>% 
  unique() %>% 
  ggplot(aes(x = conc, y = v0, color = row, shape = as.character(sensor_conc))) + 
  geom_point()
#write_tsv(processed.data, file.path(output, "fit_data.txt"))
plot_fits(processed.data, output = output)




dataset %>% 
  filter(Time < 750 & sample == "PE1_WT") %>% 
  ggplot(aes(Time, product_conc, color = conc)) + 
  geom_point() + facet_wrap(~sensor_conc) +
  scale_color_viridis()
