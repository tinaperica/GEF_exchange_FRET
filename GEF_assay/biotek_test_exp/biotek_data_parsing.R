library(lubridate) #### used the hms function for converting hh:mm:ss to seconds
library(tidyverse)
library(stringr)
data <- data.frame()
####### a function to read in all the biotek files, gather the data and reformat time
read_and_gather <- function(file) {
  data_in <-  read_delim(file, delim = "\t", col_names = T, skip = 48)
  #data_in <- select(data_in, -one_of("Temp"))
  data_in <- select(data_in, -Temp)
  data_gathered <- as.tibble(gather(data_in, key = well, value = fluorescence, -Time))
  data_gathered$Time <- as.numeric(hms(data_gathered$Time))   #### hms gives a warning when parsing 00:00:00
  relevant.index <- filter(index, data_file == file)
  dataset <- inner_join(data_gathered, relevant.index, by = "well")
  dataset <- select(dataset, -data_file)
  return(dataset)
}
# combine all the files into one tibble
(files <- dir("GEF_assay/biotek_test_exp/data/good_data", pattern = ".txt", full.names = T))
outfile <- "GEF_assay/biotek_test_exp/data/good_data_parsed.txt"
### load the index file (has conditions per well)
(index <- as.tibble(read_delim("GEF_assay/biotek_test_exp/good_data_index.txt", col_names = T, delim = "\t")))
### read in the data files, join them with the info from the index file and make them tidy
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
##### add code to remove bad data (based on a file with a list of data points to remove)
#wells_to_discard <- c()
#dataset <- filter(dataset, !(well %in% wells_to_discard) )
###############################
dataset_to_save <- mutate(dataset, "condition" = paste(date, sample, well, conc, sep = "-")) %>%
  mutate("row" = str_sub(well, 1, 1), "column" = str_sub(well, 2))
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

filter(dataset_to_plot, conc > 0.5) %>%
  arrange(conc) %>%
  ggplot(aes(x = Time, y = norm_fluorescence, color = conc)) +
    facet_wrap(~ sample) + geom_point(stroke = 0, shape = 16, size = 0.8, show.legend = T) +
  scale_colour_gradientn(colours = terrain.colors(40))
ggsave("GEF_assay/biotek_test_exp/by_sample_good_data_norm_fluorescence.pdf",  width = 15, height = 10, units = "in")


#ggplot(dataset_to_plot, aes(x = Time, y = norm_fluorescence, color = condition)) + 
 #   facet_wrap(~ sample) +
  #  geom_point(show.legend = T)
write_tsv(dataset_to_save, outfile, na = "NA", append = FALSE)
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
