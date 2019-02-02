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

file <- "GTPase_assay/Bradfords/20190201_Bradford.txt"
index <- read_tsv("GTPase_assay/Bradfords/20190201_Bradford_index.txt")
calibration_index <- index %>% filter(!is.na(conc))
data_index <- index %>% filter( is.na(conc) )

data <- read_and_gather(file)

data %>%  
  ggplot(aes(x = `450`, y = `595`)) + geom_point()
calibration <- calibration_index %>% 
  inner_join(., data, by = "well") %>% 
  mutate("conc" = as.double(conc)) %>% 
  arrange(conc)
linear_fit <- lm(conc ~ ratio_abs, data = calibration)
cal_slope <- linear_fit$coefficients[2]
cal_intercept <- linear_fit$coefficients[1]
calibration %>% ggplot(aes(x = ratio_abs, y = conc)) + 
  geom_point() + geom_abline(slope = cal_slope, intercept = cal_intercept)
gap_data <- data_index %>% 
  inner_join(., data_gathered, by = "well") %>% 
  group_by(sample) %>% summarize("mean" = mean(ratio_abs)) %>% 
  mutate("conc_ug_ml" = (cal_slope * mean + cal_intercept), "conc_uM" = conc_ug_ml / 43.25)
write_tsv(calibration, path = "GTPase_assay/Bradfords/20190201_Brad_cali.txt")
write_tsv(gap_data, path = "GTPase_assay/Bradfords/20190201_Brad_data.txt")
