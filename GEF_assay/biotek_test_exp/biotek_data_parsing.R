library(lubridate) #### used the hms function for converting hh:mm:ss to seconds
library(tidyverse)
library(stringr)
data <- data.frame()
####### a function to read in all the biotek files, gather the data and reformat time
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
    data_in <-  read_delim(file, delim = "\t", col_names = T, skip = (first_data_row - 1),
                           locale = locale(encoding = 'UTF-8'))
  }
  
  data_in <- data_in %>% select(., -2)
  data_gathered <- as.tibble(gather(data_in, key = well, value = fluorescence, -Time))
  data_gathered$Time <- as.numeric(hms(data_gathered$Time))   #### hms gives a warning when parsing 00:00:00
  relevant.index <- filter(index, data_file == file)
  dataset <- inner_join(data_gathered, relevant.index, by = "well")
  dataset <- select(dataset, -data_file)
  return(dataset)
}
# combine all the files into one tibble
(files <- dir("GEF_assay/2018_data", pattern = "GEF_FRET_assay", full.names = T))
outfile <- "GEF_assay/good_data_parsed.txt"
### load the index file (has conditions per well)
(index <- read_tsv("GEF_assay/2018_data/data_index.txt", col_names = T))
### read in the data files, join them with the info from the index file and make them tidy
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
##### add code to remove bad data (based on a file with a list of data points to remove)
#wells_to_discard <- c()
#dataset <- filter(dataset, !(well %in% wells_to_discard) )
###############################
dataset_to_save <- mutate(dataset, "condition" = paste(date, sample, well, conc, sep = "-")) %>%
  mutate("row" = str_sub(well, 1, 1), "column" = str_sub(well, 2))

write_tsv(dataset_to_save, outfile, na = "NA", append = FALSE)

dataset_to_plot <- filter(dataset_to_save, Time < 1000) %>%
  group_by(condition) %>%
  mutate("norm_fluorescence" = 1 - ((quantile(fluorescence, prob = 0.999, na.rm = T) - fluorescence) / 
                                (quantile(fluorescence, prob = 0.999, na.rm = T) - quantile(fluorescence, prob = 0.001, na.rm = T))) )

ggplot(dataset_to_plot, aes(x = Time, y = fluorescence, color = condition)) + facet_wrap(~conc) + 
  geom_point(show.legend = F)
ggsave("GEF_assay/biotek_test_exp/good_data_raw_fluorescence.pdf", width = 15, height = 10, units = "in")
ggplot(dataset_to_plot, aes(x = Time, y = norm_fluorescence, color = sample)) +
  facet_wrap(~ conc) + geom_point(stroke = 0, shape = 16, size = 0.8, show.legend = T)
ggsave("GEF_assay/biotek_test_exp/good_data_norm_fluorescence.pdf",  width = 15, height = 10, units = "in")

dataset_to_plot %>%
  arrange(as.numeric(conc)) %>%
  ggplot(aes(x = Time, y = norm_fluorescence, color = as.numeric(conc))) +
    facet_wrap(~ sample) + geom_point(stroke = 0, shape = 16, size = 0.8, show.legend = T) +
  scale_colour_gradientn(colours = terrain.colors(40))
ggsave("GEF_assay/biotek_test_exp/by_sample_good_data_norm_fluorescence.pdf",  width = 15, height = 10, units = "in")

plots <- list()
exp_dates <- dataset_to_plot %>% pull(date) %>% unique()
for (i in seq_along(exp_dates)) {
  plots[[i]] <- dataset_to_plot %>% 
    filter(date == exp_dates[i]) %>%
    ggplot(aes(x = Time, y = norm_fluorescence, color = as.numeric(conc))) + 
    facet_wrap(~ sample) + geom_point(stroke = 0, shape = 16, size = 1, show.legend = T) +
    ggtitle(exp_dates[i])
}
pdf("GEF_assay/biotek_test_exp/by_date_good_data_norm_fluorescence.pdf", width = 10, height = 5)
print(plots)
dev.off()


plots <- list()
exp_dates <- dataset_to_plot %>% 
  pull(date) %>% unique()
samples <- dataset_to_plot %>% 
  pull(sample) %>% unique()
plot_count <- 1
for (i in seq_along(exp_dates)) {
  for (j in seq_along(samples)) {
    temp <- dataset_to_plot %>% 
    filter(date == exp_dates[i], sample == samples[j])
    if (nrow(temp) > 0) {
      concentrations <- temp %>% arrange(conc) %>% pull(conc)
      temp <- temp %>% mutate("conc" = factor(conc, concentrations))
      plots[[plot_count]] <- ggplot(temp, aes(x = Time, y = norm_fluorescence, color = condition)) + 
        geom_point(show.legend = T) +
        ggtitle(exp_dates[i])
      plot_count <- plot_count + 1
    }
  }
}
pdf("GEF_assay/biotek_test_exp/by_date_and_sample_good_data_norm_fluorescence.pdf", width = 15, height = 9)
print(plots)
dev.off()

# dataset_to_save <-  select(dataset_to_save, date, Time, row, column, sample, conc, GEF_conc, fluorescence, cutoff_time)
# #write.table(dataset_to_save, file = "GEF_assay/data/20180131_biotek_data.txt", quote = F, row.names = F, sep = "\t", col.names = F)
# #### compare the reproducibility to the Tecan experiment
# (tecan_data <- as.tibble(read.delim("GEF_assay/data/GEF_FRET_kinetics.txt", head = F)))
# names(tecan_data) <- expression( date, time, row, column, sample, conc, GEF_conc, fluorescence, cutoff_time)
# tecan_data %>%
#   filter(time < 750 & (date == "20171201" | date == "20171202")) %>%
#   mutate("well" = paste0(row, column)) %>%
#   mutate("condition" = paste(well, conc, sep = "-")) %>%
#   group_by(condition) %>%
#   mutate("norm_fluorescence" = fluorescence/max(fluorescence)) %>%
#   ggplot(aes(x = time, y = norm_fluorescence, color = well)) + facet_wrap( ~conc, nrow = 4) + geom_point(show.legend = T)
