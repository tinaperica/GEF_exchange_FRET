#### import libraries
library(tidyverse)
library(viridis)
library(lubridate)
library(ggrepel)

calibration <- tibble("sensor_conc" = c(10, 20), 
                      "cal_slope" = c(0.0002841, 0.0003713), 
                      "cal_intercept" = c(-0.0014423, 0.1462813))
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
    filter(! is.na(fluorescence))
  # buffer_blanks <- data %>% filter(conc == 0) %>% 
  #   select(sensor_conc, Time, fluorescence) %>% 
  #   arrange(Time) %>% 
  #   group_by(Time) %>% 
  #   mutate("fluorescence" = mean(fluorescence)) %>% 
  #   ungroup() %>% unique()
  no_GAP_blank <- data %>% filter(GAP_conc == 0) %>% 
    select(sample, date, conc, sensor_conc, Time, fluorescence) %>% 
    filter(Time < 100) %>% 
    group_by(sample, date, conc, sensor_conc) %>% 
    summarise("mean_no_GAP_fluorescence" = mean(fluorescence)) %>% 
    ungroup()
  data <- data %>% filter(conc != 0) %>% 
    inner_join(., no_GAP_blank, by = c("sample", "date", "conc", "sensor_conc")) %>% 
    mutate("corrected_fluorescence" = fluorescence - mean_no_GAP_fluorescence) %>% 
    inner_join(., calibration, by = "sensor_conc") %>% 
    group_by(sample) %>% 
    mutate("product_conc" = as.integer(corrected_fluorescence)*cal_slope + cal_intercept) %>% 
    filter(! is.na(fluorescence)) %>% 
    ungroup()
  return(data)
}

# combine all the files into one tibble
outfiles <- "GTPase_assay/DELA_input_files"
index <- read_tsv("GTPase_assay/2018_data/data_index.txt", col_names = T)
files <- index %>% 
  pull(data_file) %>% unique()

### read in the data files, join them with the info from the index file and make them tidy
#data_points_to_discard <- read_tsv("GEF_assay/2018_data/data_to_discard.txt")
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
noGAP_blanks <- dataset %>% filter(GAP_conc == 0)
dataset <- dataset %>% filter(GAP_conc != 0)

### export data for DELA
conditions <- dataset %>% pull(condition) %>% unique()
samples <- dataset %>% pull(sample) %>% unique()
for (con in seq_along(conditions)) {
  temp <- dataset %>% 
    filter(condition == conditions[con]) %>% 
    select(condition, Time, fluorescence)
  filename <- file.path(outfiles, str_c(conditions[con], ".txt"))
  write_tsv(temp, path = filename)
}
for (s in seq_along(samples)) {
  plots <- list()
  conds <- dataset %>% filter(sample == samples[s]) %>% pull(condition) %>% unique()
  for (con in seq_along(conds)) {
    temp <- dataset %>% filter(sample == samples[s] & condition == conds[con]) %>% 
      select(condition, Time, product_conc)
    max_product_conc <- round(max(temp$product_conc, na.rm = T), 2)
    plots[[con]] <- temp %>% ggplot(aes(x = Time, y = product_conc)) + geom_point() + 
      ggtitle(str_c(conds[con], " substrate conc = ", max_product_conc, " uM"))
  }
  pdf(str_c(outfiles, "/", samples[s], ".pdf"))
  print(plots)
  dev.off()
}

#### check how much variation there is in the starting no GAP fluorescence signal
### (this is a good estimate in how much variation in signal there is due to the error in sensor amounts + variation in free Pi contamination)
### this variation will give me an overestimated error in sensor pipetting
noGAP_blanks %>% 
  select(sample, well, date, mean_no_GAP_fluorescence) %>% unique() %>% 
  ggplot(aes(x = sample, y = mean_no_GAP_fluorescence, color = as.character(date))) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


dataset %>% group_by(sample, date, well, mean_no_GAP_fluorescence) %>% 
  summarize("max_flo" = max(fluorescence, na.rm = T)) %>% 
  ungroup() %>% 
  select(sample, date, mean_no_GAP_fluorescence, max_flo) %>% 
  gather(signal, fluorescence, -sample, -date) %>% 
  ggplot(aes(x = sample, y = fluorescence, color = as.character(date))) + 
  geom_point(size = 3, alpha = 0.4) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("fluorescence range (no GAP starting fluorescence and max fluorescence)")


final_product_conc <- dataset %>% 
  group_by(sample, well, date) %>% 
  summarize("final_product_conc" = max(product_conc, na.rm = T)) %>% 
  ungroup()


parameters <- read_tsv("GTPase_assay/DELA_analysis/parameters_summary.txt")
parameters <- parameters %>% 
  separate(col = data, into= c("date", "sample", "well", "conc", "GAP_conc", "sensor_conc"), sep = "-")

avg_parameters <- parameters %>% 
  group_by(sample) %>% 
  summarize("mean_kcat" = mean(kcat), "mean_Km" = mean(Km), 
            "kcat_sd" = sd(kcat), "Km_sd" = sd(Km), 
           "sd" =  sqrt( (kcat_sd/mean_kcat)^2 + (Km_sd/mean_Km)^2 ))

#### kcat vs Km scatterplot
avg_parameters %>% 
  ggplot(aes(x = mean_kcat, y = mean_Km, label = sample)) + 
  geom_point(aes(size = sd), alpha = 0.5) +
  geom_text_repel() +
  xlab(label = "kcat / s-1") +
  ylab(label = "Km / uM")

### kcat and Km barplots
  
  

