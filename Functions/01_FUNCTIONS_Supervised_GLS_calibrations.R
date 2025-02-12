#### SET UP GEOLOCATOR ANALYSES ####

#### Updated: 23-10-2023

#### Created by Diego Vicente-Sastre
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals, 
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)


################################################################################

# FUnction to read any light data file from different models producers of 
# light-level geolocators
# @producer : brand or producer of the model
# @geo_trip : ID of the geolocator or trip
# @file_name : name of the light file to open
# @start_date : starting date the data to analyse
# @end_Date : ending date the data to analyse
# @MTexcep : IDs of those loggers from Migrate Technology, that cannot be opened
# using the readMTlux function

read_lightGLS <- function(producer, geo_trip, file_name, start_date, end_date, MTexcep) {
  
  if (file.exists(file_name)) {
    if (Producer %in% c('BAS', 'Biotrack', 'Lotek')) {
      lu <- read.table(file_name, sep = ",", skip = 1, header = FALSE, fill = TRUE)
      lu$Date <- datetime_conversion(lu[, 2])
      lu$Light <- lu[, 4]
      lu <- lu[, c(5, 6)]
    } else if (Producer == 'Migrate Technology') {
      if(geo_trip %in% MTexcep) {
        cat(paste0('File from geo ', geoID, ' cannot be opened by readMTlux function'))
        lu<-read.csv(file_name, skip = 19, sep = '\t', 
                     col.names = c('Date', 'Light'), 
                     colClasses = c('character','numeric'))
        lu$Date <- as.POSIXct(strptime(lu$Date, "%d/%m/%Y %H:%M:%S", tz = "GMT"))
        lu<-lu%>%dplyr::filter(!is.na(Date))
      } else{lu<-readMTlux(file_name)}
    }
    
    lu <- lu %>%
      mutate(dtime = as.POSIXct(Date, tz = 'GMT'),
             lux = as.numeric(gsub('\\,', '.', Light)),
             date = as.Date(Date),
             time = strftime(Date, format = "%H:%M:%S")) %>%
      filter(date >= start_date & date <= end_date)
    
    return(lu)
    
  } else {
    message(paste0('File from geo ', geoID, ' does not exist'))
    return(NULL)
  }
}

################################################################################

# Function to plot an available map inside the function BAStag::preprocessLight
# @xlim : vector includes the limits of the x axis (longitude) of the map
# @ylim : vector includes the limits of the y axis (latitude) of the map

plotMapp <- function(xlim, ylim) {
  plot.new();sp::plot(map, xlim = xlim, ylim = ylim); plot.window(xlim, ylim) 
}

################################################################################

# Function to perform the calibration.
# The calibration process consists in correct the time of the twilight events or
# delete the erroneous ones. To get this, we check all the light curves during
# the calibration period and comparing to the previous and following ones. 
# @calib_start : starting date of the calibration period
# @calib_end : ending date of the calibration period
# @threshold : light threshold of the function preprocessLight to determine 
# 'daylight' and 'night' durations.
# @lon : longitude of the place where was calibrated the logger
# @lat : latitude of the place where was calibrated the logger

calib_period<- function (calib_start, calib_end, threshold, lon, lat) {
  PD <- d.lux %>% filter(Date > calib_start & Date < calib_end)
  if (nrow(PD) == 0 || round(as.numeric(difftime(max(PD$Date), min(PD$Date), 
                                                 units = 'days')))*2 < 3) {
    cat('No calibration for step', step, '\n')
    return(c('NA', 'NA', 'NA', 'NA'))
  } else {
    lightImage(tagdata = PD, offset = 12, zlim = c(0, 20))
    tsimageDeploymentLines(PD$Date, lon = lon, lat = lat, offset = 12, lwd = 3, 
                           col = adjustcolor('orange', alpha.f = 0.5))
    tm_step <- as.POSIXct(c(calib_start, calib_end), tz = 'GMT')
    abline(v = tm_step, lwd = 2, lty = 2, col = 'orange')
    twl <- BAStag::preprocessLight(PD, threshold = threshold, offset = 12, 
                                   lmax = 12, map = TRUE, plotMap = plotMapp)
    twl <- subset(twl, Deleted == 'FALSE')
    return(twl)
  }
}

################################################################################

# Function to save trn.files from calibration periods
# @twl : twilight events in format Twilight and Rise (TRUE,FALSE)
# @geoID : ID of the logger to name the file
# @sp : species name of the logger to name the file
# @year : recovery year of the logger to name the file
# @colony : colony of the species of the logger to name the file
# @step : which calibration period is (1, 2, 3 or pre-deployment and post-recovery)


save_calib <- function(twl, geoID, sp, year, colony, step) {
  trn <- twl %>%
    mutate(date_time = Twilight,
           type = ifelse(Rise == 'TRUE', 'Sunrise', 'Sunset'),
           confidence = Marker) %>%
    dplyr::select(date_time, type, confidence)
  
  if(step == 1) step <- 'predep'
  if(step == 2) step <- 'predep'
  if(step == 3) step <- 'postrec'
  if(step == 4) step <- 'postrec'
  
  file_name <- paste(geoID, sp, year, colony, step, '.trn', sep = '_')
  write.table(trn, file_name, sep = ',', col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat('Trn saved for calibration period', step, '\n')
}

################################################################################

# Function to calculate zenith angle and twilith error distribution (alpha and sd)
# The function try to use the function thresholdCalibration and return 
# 'NA' values in case it is impossible to calculate the twilight error
# @twl : twilight events in format Twilight and type (1,2)
# @lon : longitude of the place where was calibrated the logger
# @lat : latitude of the place where was calibrated the logger
# @method : method applied in the thresholdCalibration function

calibration_SGAT <- function(twl, lon, lat, method) {
  calib_SGAT <- try(thresholdCalibration(twl$Twilight, twl$Rise, lon, lat, method = method))
  if (is.error(calib_SGAT)) {
    print('Error SGAT calibration')
    return(NA)
  } else {
    calib_SGAT <- thresholdCalibration(twl$Twilight, twl$Rise, lon, lat, method = method)
    return(calib_SGAT)
  }
}

################################################################################

# Function to check and replace outliers with NA values
# After all loggers are calibrated and twilight errors are calculated,
# we check possible outliers and uncoherent values
# @data : summary of twilight errors and zenith angles of each logger
# @var_name : vector of variables' names to calculate and replace outliers

check_and_replace_outliers <- function(data, var_names) {
  for (var_name in var_names) {
    boxplot_result <- boxplot(data[[var_name]])
    outliers <- boxplot_result$out
    data <- data %>%
      mutate({{var_name}} := ifelse(data[[var_name]] %in% outliers, NA, data[[var_name]])
      )
  }
  return(data)
}

################################################################################
################################################################################

# Function to join calibrations information dataset
# Joining of different datasets with differences between some columns.
# @file_path : path of the file contains the dataset
# @producer : producer of the loggers of the dataset
# @colony : colony where the loggers were deployed 

calib_unificated <- function(file_path, producer, colony) {
  tryCatch({
    df <- read_csv(file_path, show_col_types = FALSE)
    
    if (basename(file_path) == 'Manual_calibrations_SGAT_Sara2.csv') {
      MT_geo <- GLS.info %>% filter(Producer == 'Migrate Technology')
      df <- df %>%
        add_column(
          Producer = ifelse(df$geo_alfanum %in% MT_geo$geo_alfanum, 'Migrate Technology', 'Biotrack')
        ) %>%
        dplyr::select(-producer)%>%
        rename(Colony = colony)%>%
        relocate(Producer, .after = 'geo_alfanum')
      rm(MT_geo)
    } else if (basename(file_path) == 'Manual_calibrations_SGAT_CALEDW_updated.csv') {
      df<-df%>%
        rename(Colony = colony) %>%
        mutate(Colony = case_when(
          Colony == 'Bravadriftadj.lux' ~ 'Brava',
          Colony == 'MClaradriftadj.lux' ~ 'MClara',
          Colony == 'Raso  0.lig' ~ 'Raso',
          Colony == 'CVelho  0.lig' ~ 'CVelho',
          TRUE ~ as.character(Colony))
        )
      df <- df %>% rename(Producer = producer)
    } else if (basename(file_path) == 'Manual_calibrations_SGAT_Pablo.csv'){
      df<- df%>%
        rename(Colony = colony)%>%
        mutate(Colony = 'Veneguera')%>%
        add_column(Producer = 'Biotrack')%>%
        relocate(Producer, .after = 'geo_alfanum')%>%
        distinct()
    } else if (basename(file_path) == 'Manual_calibrations_SGAT_Biotrack_Sergi.csv'){
      df<- df%>%
        rename(Colony = colony)%>%
        mutate(Colony = 'Veneguera')%>%
        rename(Producer = producer)%>%
        relocate(Producer, .after = 'geo_alfanum')%>%
        distinct()
      
    } else if (basename(file_path) == 'Manual_calibrations_SGAT_Sergi.csv'){
      df<- df%>%
        rename(Colony = colony)%>%
        mutate(Colony = 'Veneguera')%>%
        add_column(Producer = 'Migrate Technology')%>%
        relocate(Producer, .after = 'geo_alfanum')%>%
        distinct()
    } else if (basename(file_path) == 'Manual_calibrations_SGAT_CALEDW_RdeJunco.csv'){
      df<- df%>%
        rename(Producer = producer, Colony = colony)%>%
        distinct()
    } else {
      df <- df %>%
        rename(Colony = colony)
    }
    
    return(df)
  }, error = function(e) {
    cat("Error reading or processing file:", file_path, "\n")
    return(NULL)
  })
}