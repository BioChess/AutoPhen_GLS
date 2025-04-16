#### SET UP GEOLOCATOR ANALYSES ####

#### Updated: 23-10-2023

#### Created by Diego Vicente-Sastre
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals, 
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)

################################################################################




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
