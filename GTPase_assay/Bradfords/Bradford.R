library(tidyverse)
### Bradford
## READ and FORMAT the raw ASCII data from Synergy H1
####### a function to read in all the biotek files, gather the data and reformat the time columns
read_and_gather <- function(file) {
  print(file)
  raw_in <-  read_lines(file)
  two_rows_before_data_row <- grep(raw_in, pattern = "Results")
  n_rows_to_read <- 16
  data_in <-  read_tsv(file, col_names = F, skip = (two_rows_before_data_row + 1),
                         n_max = n_rows_to_read, locale = locale(encoding = 'windows-1250'))
  names(data_in) <- c("row", as.character(seq(1, 12, 1)), "nM") 
  data_in <- data_in %>% 
    mutate("row" = c('A', 'A', 'B', 'B', 'C', 'C', 'D', 'D', 'E', 'E', 'F', 'F', 'G', 'G', 'H', 'H'))
  data_gathered <- data_in %>% 
    gather(., key = col, value = abs, -nM, -row) %>% 
    filter( ! is.na(abs)) %>% 
    mutate("well" = str_c(row, col, sep = "")) %>% 
    select(-row, -col) %>% 
    spread(key = nM, value = abs) %>% 
    mutate("ratio_abs" = as.double(`595`)/as.double(`450`))
  return(data_gathered)
}

file <- "GTPase_assay/Bradfords/20190220_Bradford.txt"
index <- read_tsv("GTPase_assay/Bradfords/20190220_Bradford_index.txt")
calibration_index <- index %>% filter(!is.na(conc))
data_index <- index %>% filter( is.na(conc) )

data <- read_and_gather(file)

data %>%  
  ggplot(aes(x = `450`, y = `595`)) + geom_point()
calibration <- calibration_index %>% 
  inner_join(., data, by = "well") %>% 
  mutate("conc" = as.double(conc)) %>% 
  arrange(conc) %>% 
  mutate("row" = str_sub(well, 1, 1), "column" = str_sub(well, 2))
low_conc <- c("A", "B", "C", "D", "E", "F", "G", "buffer")
high_conc <- c("a", "b", "c", "d", "e", "f", "g", "g_half", "h", "i", "buff")

### enzyme 
calibration_low <- calibration %>% filter(sample %in% low_conc)
linear_fit <- lm(conc ~ ratio_abs, data = calibration_low)
cal_slope <- linear_fit$coefficients[2]
m_var <- (coef(summary(linear_fit))[2, 2])^2  # slope variance
cal_intercept <- linear_fit$coefficients[1]
b_var <- (coef(summary(linear_fit))[1, 2])^2 # intercept variance
calibration_low %>% 
  ggplot(aes(x = ratio_abs, y = conc, color = row)) + 
  geom_point() + geom_abline(slope = cal_slope, intercept = cal_intercept)

gap_data <- data_index %>% 
  inner_join(., data, by = "well") %>% 
  group_by(sample) %>% 
  summarize("mean" = mean(ratio_abs), "x_var" = var(ratio_abs)) %>% 
  mutate( "ug_ml_conc_error" = sqrt( (m_var * x_var + cal_slope^2 * x_var + mean^2 * m_var) + b_var ) ) %>% 
  mutate("conc_ug_ml" = (cal_slope * mean + cal_intercept)) %>% 
  mutate("percent_error" = ug_ml_conc_error / conc_ug_ml * 100) %>% 
  mutate("conc_uM" = conc_ug_ml / 43.25, "um_error" = percent_error/100 * conc_uM)

write_tsv(calibration, path = "GTPase_assay/Bradfords/20190214_Brad_cali.txt")
write_tsv(gap_data, path = "GTPase_assay/Bradfords/20190214_Brad_data.txt")


### enzyme 
calibration_high <- calibration %>% filter(sample %in% high_conc)
linear_fit <- lm(conc ~ ratio_abs, data = calibration_high)
cal_slope <- linear_fit$coefficients[2]
m_var <- (coef(summary(linear_fit))[2, 2])^2  # slope variance
cal_intercept <- linear_fit$coefficients[1]
b_var <- (coef(summary(linear_fit))[1, 2])^2 # intercept variance
calibration_high %>% 
  ggplot(aes(x = ratio_abs, y = conc, color = row)) + 
  geom_point() + geom_abline(slope = cal_slope, intercept = cal_intercept)

ran_data <- data_index %>% 
  inner_join(., data, by = "well") %>% 
  group_by(sample) %>% 
  summarize("mean" = mean(ratio_abs), "x_var" = var(ratio_abs)) %>% 
  mutate( "ug_ml_conc_error" = sqrt( (m_var * x_var + cal_slope^2 * x_var + mean^2 * m_var) + b_var ) ) %>% 
  mutate("conc_ug_ml" = (cal_slope * mean + cal_intercept)) %>% 
  mutate("percent_error" = ug_ml_conc_error / conc_ug_ml * 100) %>% 
  mutate("conc_uM" = conc_ug_ml / 24.8, "um_error" = percent_error/100 * conc_uM)



