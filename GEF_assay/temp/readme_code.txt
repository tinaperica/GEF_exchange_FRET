get_region_to_fit <- function (data) {
  plateau_fluorescence <- quantile(data$fluorescence, prob = .7, na.rm = T)
  plateau_times <- sort(data$time[ data$fluorescence < plateau_fluorescence ])
  cutoff_time <- plateau_times[1]
  relevant_data <- data[data$time < cutoff_time,]
  relevant_data <- relevant_data[order(relevant_data$time),]
  return(relevant_data)
}
get_linear_region <- function(data) {
  data <- cbind(data, "conc_diff" = max(data$substrate_conc) - data$substrate_conc )
  max_diff <- max(data$conc_diff)
  reaction_percentages <- seq(0.1, 0.5, 0.05)
  points_check <- TRUE
  for (i in reaction_percentages) {
    if (points_check) {
      linear_times <- sort(data$time[data$conc_diff < (i * max_diff)])
      if ( length (linear_times) > 15) {
        points_check <- FALSE
      }
    }
  }
  cutoff_time <- linear_times[length(linear_times)]
  linear_data <- data[data$time < cutoff_time,]
  linear_data <- linear_data[order(linear_data$time),]
  return(linear_data)
}