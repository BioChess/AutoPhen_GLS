#### SET UP GEOLOCATOR ANALYSES ####

#### Updated: 23-10-2023

#### Created by Diego Vicente-Sastre
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals, 
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)

################################################################################

# Function to select the zenith angle
# The function selects the zenith angle determined during the calibration 
# analysis. If the geolocator was not calibrated or was not in the dataset,
# the zenith angle and the twilight error parameters (alpha mean and alpha sd) 
# are select from the summary dataset.
# We select the mean for this producer, colony and year
# @info : dataframe including the information about the logger
# @calib.summary : dataframe including the summary of zenith angle and 
# twilight error parameters 
# @calibrations : dataframe including the zenith angle and 
# twilight error parameters 

select_sun_angle <- function(producer, colony, year, calib.summary, calibrations, param = T) {

  if (info$Calib.param == 'No') {
    cat(paste0('The logger has no calibrations \n'))
    calib <- calib.summary %>%
      filter(Producer == producer, Colony == colony, Year == year)
  } else {
    cat(paste0('The logger has calibrations \n'))
    calib <- calibrations %>%
      filter(geo_alfanum == geoID)
  }
  
  if (nrow(calib) > 0) {
    zenith_mean <- ifelse(info$Calib.param == 'No', calib$zenith_mean, calib$zenith)
    zenith0_mean <- ifelse(info$Calib.param == 'No', calib$zenith0_mean, calib$zenith0)
    alpha.mean <- ifelse(info$Calib.param == 'No', calib$alpha.mean_mean, calib$alpha.mean)
    alpha.sd <- ifelse(info$Calib.param == 'No', calib$alpha.sd_mean, calib$alpha.sd)
    alpha<- c(alpha.mean, alpha.sd)
    sun_angle <- 90 - zenith_mean
    calib.param <- list(zenith_mean = zenith_mean, 
                        zenith0_mean = zenith0_mean, 
                        alpha = alpha)
  } else {
    message('No calibration data found. Sun_angle set to NA \n')
    sun_angle <- NA
  }
  
  if(param == F){return(sun_angle)}
  if(param == T){return(calib.param)}
  
}


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

read_lightGLS <- function(Producer, geo_trip, file_name, date_deployed, date_retrieved, MTexcep) {
  
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
      filter(date >= date_deployed & date <= date_retrieved)
    
    return(lu)
    
  } else {
    message(paste0('File from geo ', geoID, ' does not exist'))
    return(NULL)
  }
}

################################################################################
################################################################################

# Function to load and filter resuscitated files
loadAndFilterResuscitated <- function(file_paths) {
  
  res <- lapply(file_paths, function(path) {
    read_excel(path) %>%
      dplyr::filter(Estado != 'OK' & !is.na(Estado)) %>%
      dplyr::select(geo_trip, `Archivo trn`, Estado)
  }) %>%
    bind_rows() %>%
    dplyr::filter(Estado != 'OK' & !is.na(Estado)) %>%
    dplyr::arrange(geo_trip)
  
  return(res)
}
