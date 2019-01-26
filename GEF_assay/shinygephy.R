#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(lubridate)
library(minpack.lm)
remove(list = ls())

# create output directory
today <- gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
output <- str_c(today, '_output', '/')
dir.create(output, showWarnings = FALSE)

# load data
load('../cjm_data.Rda')
dataset <- as.tibble(dataset) %>% filter(k_est != 0)

# create output directory and files for each sample if not already made (won't overwrite curves we've already fit)
for (samp in unique(unlist(dataset$sample))) {
  
  dir.create(paste0(output,samp,'/'), showWarnings = FALSE)
  
  if(!file.exists(paste0(output, samp, '/', samp, '.csv'))) {
    data.to.write <- dataset %>%
      filter(sample == samp) %>%
      filter(Time == 0) %>%
      select(-Time, -fluorescence)
    
    data_columns <- c('max_flo', 'f_plateau', 'span1', 'span2', 'k', 'k_background','f_mid', 'f0', 'vf0')
    data.to.write[,data_columns] <- 0.0
    
    write_delim(data.to.write, path = paste0(output, samp, '/', samp, '.csv'), delim = '\t', na = 'NA')
    }
}

run_nls <- function(data, start, upper, lower, deadtime, debug = FALSE) {
  
  c0 <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  k_est <- data$k_est[1]
  #### debug mode: if debug is TRUE, print out condition being fit 
  if (debug) { print(unique(data$condition)) }

  out <- nlsLM(observed ~ span1 * exp(-k * (Time+deadtime)) + span2 * (exp(-k_background * (Time+deadtime))) + f_plateau,
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
  k_background <- coef(out)[5]

  data$predicted <- span1 * exp(-k * (data$Time+deadtime)) + span2 * (exp(-k_background * (data$Time+deadtime))) + f_plateau
  data$exchange <- span1 * exp(-k * (data$Time+deadtime))  + f_plateau + span2
  data$background <-  span2 * (exp(-k_background * (data$Time+deadtime))) + f_plateau

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

save_and_print <- function(d) {
  
  plot <- d %>% 
    ggplot(aes(x=Time, y=fluorescence)) + geom_point(color = 'black') +
      geom_line(mapping = aes(x = Time, y = predicted), color = 'red') +
      geom_line(mapping = aes(x = Time, y = exchange), color = 'green') +
      geom_line(mapping = aes(x = Time, y = background), color = 'blue')

  cond <- d$condition[1]
  samp <- d$sample[1]

  pdf(paste0(output, samp, '/', paste(cond, 'fit.pdf', sep = '_')))
  print(plot)
  dev.off()
  
  data.to.write <- read_delim(paste0(output, samp, '/', samp, '.csv'), delim = '\t')
  row.to.edit <- which(data.to.write$condition == cond)
  row.to.save <- d %>%
    select(conc, cutoff_time, max_flo, f_plateau, span1, span2, k, k_background, f_mid, f0, vf0) %>%
    head(1)
  
  data.to.write[row.to.edit,colnames(row.to.save)] <- row.to.save
  
  write_delim(data.to.write, path = paste0(output, samp, '/', samp, '.csv'), delim = '\t', na = 'NA')
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   titlePanel("GEF Analysis Fitting"),
   sidebarLayout(
      sidebarPanel(
        
        selectInput("sample", "Variant", choices = c('PE1_WT'), selectize = TRUE),
        selectInput("condition", "Curve", choices = c('20181003-PE1_WT-A5-5'), selectize = TRUE),
        sliderInput("cutoff", "Cutoff", min = 0, max = 10000, value = 5000, step = 5),
        
        sliderInput("fp_c", "fp constraints", min = 0, max = 10000, value = c(0,10000)),
        sliderInput("fp_s", "fp starting", min = 0, max = 10000, value = 3491),
        sliderInput("span1_c", "span1 constraints", min = 0, max = 10000, value = c(0, 10000)),
        sliderInput("span2_c", "span2 constraints", min = 1, max = 10000, value = c(0, 10000)),
        sliderInput("k_c", "k constraints (log scale)", min = -7, max = 1, value = c(-4,-1), step = 0.1),
        sliderInput("k_background_c", "k_background constraints (log scale)", min = -10, max = -3, value = c(-10,-4), step = 0.1)
        
      ),
      mainPanel(
        actionButton("save", "Save values and plot"),
        tableOutput('df'),
        plotOutput('plot1'),
        plotOutput('plot2')
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observe({
    available_samples <- unique(unlist(dataset$sample))
    updateSelectInput(session, 'sample', choices = available_samples, selected = NULL)
  })
  
  observe({
    d_sample <- filter(dataset, sample == input$sample)
    available_conditions <- unique(unlist(d_sample$condition))
    updateSelectInput(session, 'condition', choices = available_conditions, selected = NULL)
  })
  
  observe({
    d <- filter(dataset, sample == input$sample, condition == input$condition)
    maxf <- max(d$fluorescence, na.rm = T)
    minf <- min(d$fluorescence, na.rm = T)
    
    updateSliderInput(session, "cutoff", min = 0, max = 1.5*max(d$Time, na.rm=T), value = 1.5*max(d$Time, na.rm=T), step = 5)
    updateSliderInput(session, "fp_c", min = 0, max = maxf, value = c(0,maxf))
    updateSliderInput(session, "fp_s", min = 0, max = maxf, value = maxf/2)
    updateSliderInput(session, "span1_c", min = 0, max = maxf, value = c(0, maxf))
    updateSliderInput(session, "span2_c", min = 0, max = maxf, value = c(0, maxf))
    updateSliderInput(session, "k_c", min = -7, max = 1, value = c(-4, -1), step = 0.1)
    updateSliderInput(session, "k_background_c", min = -10, max = -3, value = c(-10, -4), step = 0.1)
  })

  fit_curve <- reactive({
    
    lower <- as.double(c(input$fp_c[1],
                          input$span1_c[1],
                          input$span2_c[1],
                          10^(input$k_c[1]),
                          10^(input$k_background_c[1])))
    
    upper <- as.double(c(input$fp_c[2],
                          input$span1_c[2],
                          input$span2_c[2],
                          10^(input$k_c[2]),
                          10^(input$k_background_c[2])))
    
    start <- list(f_plateau = input$fp_s,
                  span1 = 0.5*(upper[2]+lower[2]),
                  span2 = 0.5*(upper[3]+lower[3]),
                  k = 0.5*(upper[4]+lower[4]),
                  k_background = 0.5*(upper[5]+lower[5]))

    dataset %>%
      filter(sample == input$sample) %>%
      filter(condition == input$condition) %>%
      filter(Time < input$cutoff) %>%
      mutate(cutoff_time = input$cutoff) %>%
      group_by(sample, condition) %>%
      mutate(observed = fluorescence) %>%
      do(run_nls(., start = start, upper = upper, lower = lower, deadtime = 0, debug = F)) %>%  # fit curve
      ungroup()
  })
  
  observeEvent(input$save, {
    save_and_print(fit_curve())
  })
  
  output$plot1 <- renderPlot({
    fit_curve() %>%
      ggplot(aes(x=Time, y=fluorescence)) + geom_point(color = 'black') +
        geom_line(mapping = aes(x = Time, y = predicted), color = 'red') +
        geom_line(mapping = aes(x = Time, y = exchange), color = 'green') +
        geom_line(mapping = aes(x = Time, y = background), color = 'blue')
  })
  
  output$plot2 <- renderPlot({
    fit_curve() %>%
      ggplot(aes(x=Time, y=fluorescence)) + geom_point(color = 'black') +
        geom_line(mapping = aes(x = Time, y = predicted), color = 'red')
  })
  
  output$df <- renderTable({
    fit_curve() %>%
      select(conc, cutoff_time, max_flo, f_plateau, span1, span2, k, k_background, f_mid, f0, vf0) %>%
      head(1)
  }, digits=4)
  
}


# Run the application 
shinyApp(ui = ui, server = server)

