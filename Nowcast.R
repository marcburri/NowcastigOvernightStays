#--------------------------------------------------------
####                    3. OS                       ####
#--------------------------------------------------------
# This script estimates forecasting models for the number
# of Overnight stays in 24 EU countries. This approach 
# reproduces entry 4, which won the second price in the 
# European Nowcasting Accuracy Challenge, as well as the 
# first price in the European Nowcasting Reproducibility 
# Challenge
# - Models: STLM
#--------------------------------------------------------

#### Preliminaries ####
#--------------------------------------------------------
# To be defined
nowcastDate <- as.Date("2024-03-01") # Set the date for nowcasting to March 1, 2024
testSet <- F # Boolean to choose if the test set should be forecasted (F = no, T = yes)
excludeCovid <- T # Boolean to decide if observations during Covid lockdowns should be excluded (T = yes)
selectedModels <- c( "stlm") # Specify the forecasting model to be used; here it's STLM

#--------------------------------------------------------
# Load necessary R packages for data manipulation, time series analysis, and visualization
library(eurostat) # For fetching Eurostat data
library(dplyr) # For data manipulation
library(tsbox) # For time series conversion
library(jsonlite) # For reading and writing JSON files
library(tidymodels) # For modeling framework
library(modeltime) # For time series forecasting models
library(timetk) # For time series manipulation and visualization
library(lubridate) # For date-time manipulation
library(imputeTS) # For missing data imputation in time series

# Define a tibble with acronyms of 24 EU countries
countries <- tibble(
  acronym = c(
    "AT", # Austria
    "BE", # Belgium
    "BG", # Bulgaria
    "CY", # Cyprus
    "CZ", # Czech Republic
    "DE", # Germany
    "DK", # Denmark
    "EE", # Estonia
    "EL", # Greece
    "ES", # Spain
    "FI", # Finland
    "FR", # France
    "HR", # Croatia
    "HU", # Hungary
    "IE", # Ireland
    "IT", # Italy
    "LT", # Lithuania
    "LU", # Luxembourg
    "LV", # Latvia
    "MT", # Malta
    "NL", # Netherlands
    "PL", # Poland
    "PT", # Portugal
    "RO", # Romania
    "SI", # Slovenia
    "SK", # Slovakia
    "SE"  # Sweden	
  )
)

#--------------------------------------------------------
# Create a copy of the json file to write results
results <- read_json(path = "Results/starting_kit_os/point_estimates.json") # Read the initial JSON file with results
write_json(results, pretty = T, paste0("Results/point_estimates_OS_m", month(nowcastDate),".json")) # Write the results to a new JSON file with a name based on the nowcast month

# Fetch data on Overnight Stays (OS) from the Eurostat database
OS <- get_eurostat(id="tour_occ_nim")

# Prepare the data: rename columns, filter relevant data, and select necessary columns
OS <- OS %>%
  rename(time = TIME_PERIOD) %>%
  filter(
    geo %in% countries$acronym, # Filter for the selected countries
    unit=="NR", # Unit of measurement
    nace_r2=="I551-I553", # NACE classification for accommodation
    c_resid == "TOTAL" # Total residents and non-residents
  ) %>%
  rename(
    id = geo, # Rename 'geo' column to 'id'
    OS = values # Rename 'values' column to 'OS'
  ) %>%
  dplyr::select(time, id, OS) %>% # Select only time, id, and OS columns
  arrange(time, decreasing = F) # Arrange the data by time in ascending order

#### Modelling & Forecasting ####
#--------------------------------------------------------
for (cntry in countries$acronym) { # Iterate over each country in the 'countries' tibble
  # Check if the country exists in the json file, skip to next iteration if not
  if (is.na(names(results$entry_1[cntry]))) {
    next
  }
  
  # Initialize an empty modeltime table for storing models
  models_tbl<- modeltime_table()
  
  # Set the start date for model estimation, varying by country
  if (cntry %in% c("LT", "MT")) {
    startDate <- as.Date("2005-01-01") # Later start date for Lithuania and Malta
  } else {
    startDate <- as.Date("2004-01-01") # General start date for other countries
  }
  
  # Prepare the dependent variable (OS data) for the current country
  dep <- OS %>% filter(
    id == cntry
  ) %>%
    ts_xts() %>%
    ts_tbl()
  
  # Determine the split date for the data
  splitDate <- as.Date(dep[nrow(dep),]$time)
  
  # Impute missing values and handle zero values in the data
  dep <- OS %>%
    filter(
      id == cntry,
      time <= splitDate
    ) %>%
    select(time, OS) %>%
    rename(value = OS) %>%
    ts_span(end=splitDate, extend = T) %>%
    mutate(
      value = na_interpolation(value, "stine"), # Impute missing values using Stine method
      value = if_else(value==0,1,value)) %>% # Replace zero values with 1
    ts_xts() %>%
    ts_tbl()
  
  # Exclude data from the COVID-19 period if specified
  if (excludeCovid) {
    depFcst <- dep %>% # Exclude covid crisis period
      filter(!between(time, as.Date("2019-03-01"), as.Date("2022-02-28")) ) %>%
      mutate(time = if_else(time <= as.Date("2019-03-01"), time %m+% years(3), time)) 
  } else {
    depFcst <- dep
  }
  
  # Check if forecasting beyond available data is needed and prepare for it
  if (interval(splitDate, nowcastDate) %/% months(1) > 0) {
    recurs <- T
    recursPeriods <- interval(splitDate, nowcastDate) %/% months(1)
    lastVal <- dep[nrow(dep),]$value
    dep <- ts_c(dep, ts_forecast(depFcst, recursPeriods)) %>% select(time, value)
  } 
  
  # Log-transform the dependent variable for modeling
  dta_tbl <- ts_tbl(dep) %>%
    mutate(value = log(value)) %>%
    rename(date=time)
  
  # Optionally exclude COVID-19 period data for the forecasting dataset
  if (excludeCovid) {
    dta_tbl <- dta_tbl %>%
      filter(!between(date, as.Date("2019-03-01"), as.Date("2022-02-28")) ) %>%
      mutate(date = if_else(date <= as.Date("2019-03-01"), date %m+% years(3), date)) 
  }
  
  # Prepare the forecast dataset for future dates
  dta_tbl_fcst <- dta_tbl %>%
    filter(date > splitDate)
  
  # Prepare the training dataset for model fitting
  dta_tbl_split <- dta_tbl %>%
    filter(date <= splitDate)
  
  # Create time series splits for model training and testing
  splits <- dta_tbl_split %>%
    time_series_split(date, 
                      assess = "5 months", 
                      cumulative = T)
  
  # Define a recipe for the model specifying the formula and data
  recipe_spec <- recipe(value ~ ., data = training(splits))
  
  # If STLM model is selected, fit the model to the training data
  if ("stlm" %in% c(selectedModels)) {
    wflw_fit_stlm <- workflow() %>% 
      add_model(seasonal_reg(seasonal_period_1 = 12) %>%
                  set_engine("stlm_ets")) %>%
      add_recipe(recipe_spec %>%
                   step_select(all_of(c("date", "value")))
      ) %>%
      fit(training(splits))
    
    models_tbl<- add_modeltime_model(models_tbl, wflw_fit_stlm) %>%
      update_modeltime_description(nrow(models_tbl)+1, "stlm")
  }
  
  # Calibrate the models on the test set
  calibration_tbl <- models_tbl %>% 
    modeltime_calibrate(testing(splits))
  
  # If forecasting for the test set, generate forecasts and calculate accuracy metrics
  if (testSet) {
    # Generate forecasts for the test set
    forecastTest <- calibration_tbl %>%
      modeltime_forecast(
        new_data     = testing(splits),
        actual_data = dta_tbl,
        keep_data   = TRUE 
      ) 
    
    # Plot and save the forecast vs test set comparison
    p <- forecastTest %>%
      filter(.index > as.Date("2019-01-01")) %>%
      plot_modeltime_forecast(
        .title = paste0("OS - ", cntry),
        .interactive = F
      )
    
    ggsave(p, file=paste0("Results/OS/forecast_test_",cntry,".pdf"), width = 8, height= 5)
    
    # Calculate and save accuracy metrics
    metricTest <- calibration_tbl %>%
      modeltime_accuracy(
        testing(splits),
        metric_set = yardstick::metric_set(mae, smape, rmse, rsq )) 
    
    print(xtable(metricTest, type = "latex"), file = paste0("Results/OS/forecast_test_",cntry,".tex"))
  }
  
  # Refit the model on the full dataset and generate final forecasts
  refit_tbl <- calibration_tbl %>%
    modeltime_refit(dta_tbl_split) %>%
    mutate(.model_desc = calibration_tbl$.model_desc) # Keep the model description
  
  # Calculate forecast for the future period
  forecast <- refit_tbl %>%
    modeltime_forecast(
      new_data    = dta_tbl_fcst,
      actual_data = dta_tbl,
      keep_data   = F 
    ) 
  
  # Plot and save the final forecast
  p <- forecast %>% 
    filter(.index > as.Date("2019-01-01")) %>%
    plot_modeltime_forecast(
      .title = paste0("OS - ", cntry),
      .interactive = F
    )
  
  ggsave(p, file=paste0("Results/OS/forecast_",cntry,"_m",month(nowcastDate),".pdf"), width = 8, height= 5)
  
  # Update the results JSON file with the new forecasts
  results <- read_json(path =  paste0("Results/point_estimates_OS_m", month(nowcastDate),".json"), auto_unbox = T)
  results$entry_4[[cntry]] <- round(exp(forecast %>% filter(.index == nowcastDate, .model_desc == selectedModels[1]) %>% pull(.value)),1)
  write_json(results, pretty = T,  paste0("Results/point_estimates_OS_m", month(nowcastDate),".json"), auto_unbox = T)
}