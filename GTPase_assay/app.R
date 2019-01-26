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
library(minpack.lm)
remove(list = ls())

# create output directory
today <- gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
output <- str_c(today, '_output', '/')
dir.create(output, showWarnings = FALSE)

# load data
dataset <- read_tsv("GAP_assay_data_parsed.txt")
#data <- dataset %>% filter(condition == "20190122-PE1_WT-C3-0.75-1-10")
#upper <- c(300, 2500, 1)
#lower <- c(-10, 1100, 0.0001)
run_fit <- function (data, fit, start, upper, lower, f0_intercept, percent_fit) {
  c0 <- data$conc[1]
  GAP_conc <- data$GAP_conc[1]
  if (fit == "exp") {
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
    data$k_pval <- log10(summary(out)$coefficients[,4][3])
    ### calculate the initial rate
    data$v0 = (span * k) / (GAP_conc * 0.001) ### initial rate in fluorescence units
    data$slope <- NA
    data$percent_fit <- NA
    data <- data %>% 
      mutate("initial" = "all")
  } else if (fit == "lin") {
    f0 <- f0_intercept
    max_product <- max(data$product_conc, na.rm = T)
    min_product <- min(data$product_conc, na.rm = T)
    delta_product <- max_product - min_product
    cutoff_product <- min_product + (percent_fit * delta_product)
    initial_data <- data %>% filter(product_conc < cutoff_product)
    fit.lin <- lm(I(product_conc - f0) ~ Time - 1, data = initial_data)
    slope <- fit.lin$coefficients
    data$predicted <- slope * data$Time + f0  
    data$slope <- slope
    data$f0 <- f0
    data$f_plateau <- max_product
    data$v0 <- - slope / (GAP_conc * 0.001)
    data$percent_fit <- percent_fit
    ## label the initial points used for fitting
    data <- data %>% 
      mutate("initial" = ifelse(product_conc < cutoff_product, "initial", "all")) %>% 
      mutate("predicted" = ifelse(predicted < max_product, predicted, NA))
    data$k <- NA
    data$span <- NA
    data$chisq <- NA
    data$pearsonr <- NA
    data$k_pval <- NA

  }
  return(data)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("GAP Analysis Fitting"),
  sidebarLayout(
    sidebarPanel(
      
      selectInput("sample", "Variant", choices = c('PE1_WT'), selectize = TRUE),
      selectInput("condition", "Curve", choices = c('20190122-PE1_WT-C1-0.25-1-10'), selectize = TRUE),
      selectInput("fit", "Fitting method", choices = c('exp', 'lin'), selectize = TRUE),
      sliderInput("cutoff_time", "Cutoff time", min = 0, max = 10000, value = 5000, step = 5),
      sliderInput("f0_c", "f0 constraints", min = -1, max = 2, value = c(-1, 2), step = 0.001),
      sliderInput("f_plat_c", "f_plateau constraints", min = 0, max = 50, value = c(0, 50), step = 0.001),
      sliderInput("k_c", "k constraints (log scale)", min = -7, max = 1, value = c(-4,-1), step = 0.1),
      sliderInput("percent_fit", "percent start reaction to fit", min = 0.01, max = 0.4, value = 0.01, step = 0.01)
    ),
    mainPanel(
      tableOutput('df'),
      plotOutput('plot1')
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
    updateSelectInput(session, 'fit', choices = c('exp', 'lin'), selected = NULL)
  })
  
  observe({
    d <- filter(dataset, sample == input$sample, condition == input$condition)
    max_product <- max(d$product_conc, na.rm = T)
    min_product_f0 <- round(min(1.3 * min(d$product_conc, na.rm = T) , 0.5 * min(d$product_conc, na.rm = T) , - 0.1), 3)
    max_product_f0 <- round(max(1.3 * min(d$product_conc, na.rm = T) , 0.5 * min(d$product_conc, na.rm = T) , 0.1), 3)
    updateSliderInput(session, "cutoff_time", min = 0, max = max(d$Time, na.rm=T), value = max(d$Time, na.rm = T), step = 5)
    updateSliderInput(session, "f0_c", min = -1, max = 2, value = c(min_product_f0, max_product_f0))
    updateSliderInput(session, "f_plat_c", min = round(0.5 * max_product, 3), max = round(2 * max_product, 3), value = c(round(0.5 * max_product, 3), round(2 * max_product, 3)))
    updateSliderInput(session, "k_c", min = -7, max = 1, value = c(-4, -1), step = 0.1)
    updateSliderInput(session, "percent_fit", min = 0.01, max = 0.4, value = 0.1, step = 0.01)
  })
  
  fit_curve <- reactive({
    
    lower <- as.double(c(input$f0_c[1],
                         input$f_plat_c[1],
                         10^(input$k_c[1])
                         ))
    
    upper <- as.double(c(input$f0_c[2],
                         input$f_plat_c[2],
                         10^(input$k_c[2])
                         ))
    
    start <- list(f0 = mean(c(upper[1], lower[1])),
                  f_plat = mean(c(upper[2], lower[2])),
                  k = 0.1*(upper[3]+lower[3])
                  )
    
    dataset %>%
      filter(sample == input$sample) %>%
      filter(condition == input$condition) %>%
      filter(Time < input$cutoff_time) %>% 
      do(run_fit(., fit = input$fit, start = start, upper = upper, lower = lower, f0_intercept = 0.5*(upper[1]+lower[1]), input$percent_fit)) # fit curve or line
    
  })
  
  output$plot1 <- renderPlot({
    fit_curve() %>%
      ggplot(aes(x = Time, y = product_conc, color = initial)) + geom_point() +
      geom_line(mapping = aes(x = Time, y = predicted), color = 'green')
  })
  
  
  output$df <- renderTable({
    fit_curve() %>%
      select(chisq, pearsonr, conc, f_plateau, span, k, k_pval, slope, f0, v0) %>%
      head(1)
  }, digits=4)
  
}


# Run the application 
shinyApp(ui = ui, server = server)



# test
# samp <- "PE1_WT"
# con <- "20190122-PE1_WT-C1-0.25-1-10"
# fit <- "exp"
# cutoff_time <- 5000
# f0 <- 0.5 * (1000 - 800)
# d <- dataset %>% filter(sample == samp, condition == con)
# maxf <- max(d$product_conc, na.rm = T)
# minf <- min(d$product_conc, na.rm = T)
# f0 <- 0.5 * ( (min(1.3 * minf, 0.7 * minf, 0)) + max(1.3 * minf, 0.7 * minf, 0))
# f_plat <- 0.5 * (0.5 * maxf +  1.5 * maxf)
# k <- 0.1*(10^-4 + 10^-1)
# 
# lower <- c(-800,
#                     0.5* maxf,
#                      10^-4
# )
# 
# upper <- c(1000,
#                     1.5 * maxf,
#                      10^-1
# )
# 
# start <- list(f0 = f0,
#               f_plat = f_plat,
#               k = k
# )
# data <- d
# res <- run_fit(data = data,  fit = fit, start = start, upper = upper, lower = lower, cutoff_time = cutoff_time, f0_intercept = 0.5*(upper[1]+lower[1]), percent_fit = 0.1)
