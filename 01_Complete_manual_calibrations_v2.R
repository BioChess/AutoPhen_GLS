## SEMI-AUTOMATED DEFINITION OF GLS CALIBRATIONS ##

# Script developed by Diego Vicente-Sastre
# Based on previous script developed by Teresa Militao and RaUl Ramos

# Updated: 22 - 12- 2022

#There are modifications in get.elevation function

#### 0.1 Working directory ####
Sys.setenv(TZ = 'GMT')


#### 0.2 install packages ####
#install.packages("devtools")
#library(devtools)
#install.packages("remotes")
#library(remotes)
#install_github("SWotherspoon/SGAT")
#install_github("SWotherspoon/BAStag")
#install_github("bmcclintock/momentuHMM")
#install.packages("momentuHMM")
#install_github("SLisovski/TwGeos")
#install_github("SLisovski/GeoLight")

# Installing all the packages and load their libraries we need to run the script

pacman::p_load('devtools','remotes', # loading packages
               'readxl', # opening Excels and csv
               'tidyverse',  #manging data
               'lubridate', # dealing with dates
               'maps', 'rworldmap', 'maptools', # creating map for preprocessLight
               'GeoLight', 'BAStag', 'TwGeos','SGAT' # geolocators
               )
detach("package:purrr", unload = TRUE)

#This function is loaded from scripts of SEATRACK project
functions.dir<- paste0(wd,'/RCode/Migratory patterns/SEATRACK/functions/')
source(paste0(functions.dir,'FUNCTION_datetime_conversion.R'))
source(paste0(wd,'/RCode/Migratory patterns/00_Setup_GLS.R'))
source(paste0(wd,'/RCode/Migratory patterns/02_FUNCTIONS_Manual_calibrations.R'))

#### 1. Load GLS info  ####

GLS.files <- read_csv('RawData/GLS.files.csv', col_names = TRUE)
GLS.info <- read_csv('RawData/GLS.info.csv', col_names = TRUE)
calib.info <- read_csv('RawData/Info.calibrations.csv')

col.name <- 'Veneguera'
col.lat <- 27.84444
col.lon <- -15.78861


#### 2. Loop processtwilights and calculate calibration parameters ####
# defines the light leves above which it is day and below it is night 
# (better use 1.5 for Intigeo tags if no strong reason for other value)

threshold = 2 

xlim <- c(-20, 0)
ylim <- c(15, 40)
map <- rworldmap::getMap(resolution = "coarse")

calib.list_SGAT<-list()


# Loop to process .lux to .trn and determine calibrations parameters for each GLS
for (i in 1:nrow(GLS.files)) {
  #i=1
  setwd(paste0(wd, '/Raw Data/GLS/lux/'))
  
  # Load .lux file and extract the GLS id from the name of the file
  info<-GLS.files[i,]
  geoID <- info$geo_alfanum#Which GLS is processing
  
  lightfile <- if_else(!is.na(info$File.adjlux), info$File.adjlux, info$File.lux)
  if (!file.exists(lightfile) || geoID == GLS.files[i - 1, 'geo_alfanum']) {
    message(' Skipping geo (',i,') ', geoID, ' as there is no file or it has already been calibrated.\n')
    next
  }
  
  geoID <- info$geo_alfanum#Which GLS is processing
  year <- info$Year_rec
  cat('Processing geo (', i,') ', geoID, 'for year', year, '\n')
  
  # Select colony and determine coordinates
  sp <- info$Sp
  
  
  # Get info of the GLS and its calibration

  calib<-calib.info%>%filter(geo_alfanum == geoID)# calibration info
  
  start_date<-as.Date(info$start_date, tz = 'GMT') #date start the trip
  end_date <-as.Date(info$end_date, tz = 'GMT') #date end the trip
  
  LoggerModel<-info$Modelo#logger model
  Producer<- info$Producer #e.g "Migrate Technology","BAS","Biotrack","Lotek"

  
  lu <- read_lightGLS(producer = Producer, geo_trip, 
                      file_name = lightfile, 
                      start_date = start_date, 
                      end_date = end_date, 
                      MTexcep = MTlux_exception)
  
  # d.lux$Date <- as.POSIXct(strptime(d.lux$Date, format = '%Y-%m-%d %H:%M:%S', tz = 'GMT'))

  # Calculate zenith angle and alpha parameter from calibration periods. 
  # Our devices have been calibrated as maxim three times, so we repeat the loop 
  # three times
  
  twl_error.list <-list()
  for (step in 1:3) {
  
    # Define the calib_start and calib_end for the current step
    calib_start <- calib[[paste0("calib_start", step)]]
    calib_end <- calib[[paste0("calib_end", step)]] + 1 * 60 * 60 * 24
    n_days <- round(as.numeric(difftime(calib_end, calib_start, 
                                        units = 'days')))
    loc <- calib[[paste0("loc_calib", step)]]
    lon <- calib[[paste0("Lon", step)]]
    lat <- calib[[paste0("Lat", step)]]
    
    calib_SGAT <- NA
    
    if (!is.na(loc) && !is.na(calib_end)) {
      # Perform calibration for the current step
      twl <- calib_period(calib_start = calib_start, 
                          calib_end = calib_end, 
                          lon = lon, lat = lat,
                          threshold =threshold)
      
      # Save calibration data
      save_calib(twl, geoID, sp, year, col.name, step)
      
      # Perform SGAT calibration
      calib_SGAT <- calibration_SGAT(twl, lon = lon, lat = lat, method = 'gamma')
      
    } else {
      message('No calibration information for step',step,'\n')
      calib_SGAT <- c('NA', 'NA', 'NA', 'NA')
    }
    
    if(step == 3){cat('Geo done')}
    
    twl_error <- bind_cols(geoID = geoID, Producer = Producer, Year = year, Colony = col.name,
                           n_days = n_days,
                           zenith = calib_SGAT[1], 
                           zenith0 = calib_SGAT[2], 
                           alpha.mean = calib_SGAT[3], 
                           alpha.sd = calib_SGAT[4],
                           Loc = loc)
   
    colnames(twl_error)[5:10]<- paste0(colnames(twl_error[5:10]),'_PD', step)
    
    
    twl_error.list[[step]] <-list(twl_error)
  }
  twl_error_GLS<- bind_rows(calib_SGAT.list)%>%
    group_by(geoID, Producer, year, Colony) %>% 
    summarise_all(funs(na.omit(.)[1]))
  
  calib.list_SGAT[i]<-list(twl_error_GLS)
  
}


#SGAT calibration parameters
calib.all_SGAT <- bind_rows(calib.list_SGAT)# transform list to data.frame

str(calib.all_SGAT)# check the structure of our data
calib.all_SGAT<-calib.all_SGAT%>%
  mutate_at(c(3,5:8,10:13,15:18), as.numeric)%>% # transform to numeric, then we can operate
  rename(n_days_PR3 = n_days_PD3,
         zenith_PR3 = zenith_PD3, zenith0_PR3 = zenith0_PD3, 
         alpha.mean_PR3 = alpha.mean_PD3, 
         alpha.sd_PR3 = alpha.sd_PD3,)
  mutate(zenith = NA, zenith0 = NA, alpha.mean = NA, alpha.sd = NA)# create summary columns


#### 3. Check the values of each variable of SGAT, check if there are outliers ####

# Define variables and producers to check
variables_to_check <- c("zenith_PD1", "zenith_PD2", "zenith_PR3", 
                        "zenith0_PD1", "zenith0_PD2", "zenith0_PR3", 
                        "alpha.mean_PD1", "alpha.mean_PD2", "alpha.mean_PR3", 
                        "alpha.sd_PD1", "alpha.sd_PD2", "alpha.sd_PR3")

producers_to_check <- c("Migrate Technology", "Biotrack")

calib_list <- list()

# Loop through producers and apply the outlier detection and replacement
for (producer in producers_to_check) {
  calib_data <- calib.all_SGAT %>%
    filter(Producer == producer)
  calib_data <- check_and_replace_outliers(calib_data, variables_to_check)
  calib_list[[producer]] <- calib_data
}

# Combine data from different producers
calib.all <- bind_rows(calib_list)

# Sort the data by geo_alfanum
calib.all <- calib.all %>% arrange(geo_alfanum)

#### 4. Calculate mean values of calibration parameters ####

# Group the data by year, colony, and producer
calib.all_SGAT2<- calib.all %>%
  group_by(year, Colony, Producer)%>%
  mutate(mean_zenith = ifelse(is.na(zenith_PD1) | is.na(zenith_PD2) | is.na(zenith_PR3),
                              mean(c(zenith_PD1, zenith_PD2, zenith_PR3), na.rm = TRUE),
                              zenith_PD1))%>%
  mutate(mean_zenith0 = ifelse(is.na(zenith0_PD1) | is.na(zenith0_PD2) | is.na(zenith0_PR3),
                               mean(c(zenith0_PD1, zenith0_PD2, zenith0_PR3), na.rm = TRUE),
                               zenith0_PD1))%>%
  mutate(mean_alpha.mean = ifelse(is.na(alpha.mean_PD1) | is.na(alpha.mean_PD2) | is.na(alpha.mean_PR3),
                                  mean(c(alpha.mean_PD1, alpha.mean_PD2, alpha.mean_PR3), na.rm = TRUE),
                                  alpha.mean_PD1)) %>%
  mutate(mean_alpha.sd = ifelse(is.na(alpha.sd_PD1) | is.na(alpha.sd_PD2) | is.na(alpha.sd_PR3),
                                mean(c(alpha.sd_PD1, alpha.sd_PD2, alpha.sd_PR3), na.rm = TRUE),
                                alpha.sd_PD1))%>%
  filter(!duplicated(geo_alfanum))



par(mfrow= c(3,2))
boxplot(calib.all_SGAT2[,20:21])
boxplot(calib.all_SGAT2[,22:23])
hist(calib.all_SGAT2$zenith)
hist(calib.all_SGAT2$zenith0)
hist(calib.all_SGAT2$alpha.mean)
hist(calib.all_SGAT2$alpha.sd)

calib_SGAT_summary <- calib.all_SGAT2 %>%
  # filter(
  #   mean_zenith != 'NaN' &
  #   mean_zenith0 != 'NaN' &
  #   mean_alpha.mean != 'NaN' &
  #   mean_alpha.sd != 'NaN'
  # ) %>%
  group_by(Colony, year, Producer) %>%
  summarize(
    n_geos = n_distinct(geo_alfanum),
    zenith_mean = mean(mean_zenith, na.rm = TRUE),
    zenith0_mean = mean(mean_zenith0, na.rm = TRUE),
    alpha.mean_mean = mean(mean_alpha.mean, na.rm = TRUE),
    alpha.sd_mean = mean(mean_alpha.sd, na.rm = TRUE)
  ) %>%
  arrange(year) %>%
  distinct()



#### 5. Save final results ####
write_csv(calib_SGAT_summary, 'Raw Data/Summary_SGAT_calibrations_ICELAND.csv')
write_csv(calib.all_SGAT, 'Raw Data/Manual_calibrations_SGAT_ICELAND.csv')
write_csv(calib.all_SGAT2, 'Raw Data/Manual_calibrations_SGAT2_ICELAND.csv')











