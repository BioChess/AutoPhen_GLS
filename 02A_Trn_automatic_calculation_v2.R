#### TRN AUTOMATIC CALCULATION ####

# Updated: 05-03-2023

#### Created by Diego Vicente-Sastre
#### Based on previous script developed by SEATRACK project
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals, 
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)

################################################################################

############################################################################################################

#### 0.1 Working directory ####
Sys.setenv(TZ = 'GMT') 
### !!! IMPORTANT: IF NOT SET LIKE THIS, MAKE PROBLEMS TO CREATE DATE_TIME FROM 
# PASTING DATE plus TIME

place <- 'PCUB'  # pc in the UB
#place <- "home"  # laptop at home


if (place == 'PCUB') WD <- 'D:/Dropbox/'
if (place == 'home') WD <- 'F:/Dropbox/'  # pc in the UB
setwd(paste0(WD, 'Diego Vicente/TESIS/01_GLS_analysis/'))

wd<-getwd()

functions.dir<- paste0(wd,'/RCode/Migratory patterns/SEATRACK/functions/')

#### 0.2 install packages ####

# Installing all the packages and load their libraries we need to run the script

pacman::p_load('devtools','remotes', # loading packages
               'readxl', # opening Excels and csv
               'tidyverse',  #manging data
               'lubridate', # dealing with dates
               'maps', 'rworldmap', 'maptools', # creating map for preprocessLight
               'GeoLight', 'BAStag', 'TwGeos','SGAT', # geolocators
               'suncalc', # knowing twilights depending on the place
               'tictoc')


#These function are loaded from scripts of SEATRACK project
source(paste0(functions.dir,"FUNCTION_datetime_conversion.R"))
source(paste0(functions.dir,"FUNCTION_twilight_cleanup.R"))
source(paste0(functions.dir,"FUNCTION_twilight_cleanup_v2.R"))
source(paste0(functions.dir,"FUNCTION_move_twilights.R"))
source(paste0(functions.dir,"FUNCTION_double_smoothing.R"))
source(paste0(functions.dir,"FUNCTION_assign_equinox_periods.R"))

source(paste0(wd,'/RCode/Migratory patterns/00_Setup_GLS.R'))
source(paste0(wd,'/RCode/Migratory patterns/03_FUNCTIONS_Trn_automatic_calculation.R'))
############################################################################################################

#### 1. Load GLS info  ####

GLS.files <- read_csv('Raw Data/GLS.files.csv', col_names = TRUE)
GLS.files2 <- GLS.files%>%
  add_column(File.trn = NA, No_data = NA, Do_again = NA )
#GLS.files2[1485,'geo_trip'] = 'BH062001001'
#GLS.files2[1509,'geo_trip'] = 'BZ599001001'
GLS.info <- loadGLSInfo() 

colonies <- loadCol.Coord() # colony's coordinates

calib.info <- loadCalib.Info() # calibration periods
calibrations <- loadCalib.Angles() # calibration parameters calculated
calib.summary <- loadCalib.Summary() # calibration summary per year and colony

#Equinox vary a little between years, we selecte the precise dates for each GLS
equinox<-read_excel('Raw Data/Equinox_by_year.xlsx')
head(equinox)

# Define threshold and offset
offset = 12 # adjusts the y-axis to put night (dark shades) in the middle 
LightThreshold <- 1.5

################################################################################

# Loop to get trn files for every GLS

setwd(paste0(wd, '/Raw Data/GLS/lux/'))
trn.list1<-list()
trn.list2<-list()

no_file <- c('BD797002', 'BR670001', 'BX609001', 'CE839001','CE864001')
MTlux_exception <- c('N643003001', 'Q666002001', 'X250001001')
months_breeding=c(5,6,7,8)

tic()

for (i in 1:nrow(GLS.files2)) {
  #i=1
  # i=which(GLS.files$geo_trip == '1652005001')
  
  # 3. Extract the GLS information ####
  info<-GLS.files2[i,]
  geoID<-info$geo_alfanum #Which GLS is processing
  lightfile <- if_else(!is.na(info$File.adjlux), info$File.adjlux, info$File.lux)
  N_trip <- info$N_trips
  
  if(geo_trip == 'Q666002001') lightfile<- info$File.lux #adj.lux it is bad in this GLS
  
  # We skip the GLS that do not have light file
  if (!file.exists(lightfile) | 
      is.na(lightfile) | geoID %in% no_file | N_trip == 0){
   
    if(N_trip == 0){ #Geolocators with data only from breeding period are skipped
      message('Geo (',i,') ', geoID, 'has data only from breeding period - Skipped \n')
    } else {
      message(paste0('File from ', geoID, ' does not exist - Skipped'))
    }
    
    GLS.files2[i, c(21:23)] <- list(NA, 'No data', 'No')
    next
  } 

  geo_trip<-info$geo_trip #which GLS trip is selecting
  Res <- info$Resucitat

  # We mark the resuscitated GLS, we will do apart
  if(Res == 'Yes'){
    cat('Geo (',i,') ', geoID, 'is a resuscitated file \n')
    GLS.files2[i, c(21:23)] <- list(NA, 'Data resucitated', NA)
  }
  
  LoggerModel<-info$Modelo #logger model
  Producer<- info$Producer #company
  
  year_rec <- info$Year_rec #which year we recovered the GLS
  colony <- info$Colony #which colony is from
  
  if(colony == 'M.Clara') colony <- 'MClara'
  if(colony == 'Cala Morell') colony <- 'CalaMorell'
  if(colony == 'Curral Velho') colony <- 'CVelho'
  if(colony == 'Ilheu de Cima') colony <- 'ICima'
  cat('Processing geo (', i,') ', geo_trip, 'for year', year_rec, ' from ', colony,  '\n')

  col_lon <- as.numeric(colonies[(colonies$Colony == colony),]$Lon) #colony Longitude
  col_lat <- as.numeric(colonies[(colonies$Colony == colony),]$Lat) #colony longitude
  
  sun_angle <- select_sun_angle(producer = Producer, colony, year = year_rec,
                                calib.summary, calibrations, param = F)
  
  start_date<-as.Date(info$start_date, tz = 'GMT') #date start the trip
  end_date <-as.Date(info$end_date, tz = 'GMT') #date end the trip

  
  # 3.1. Reading light data ####
  
  lu <- read_lightGLS(Producer, geo_trip, 
                      file_name = lightfile, 
                      date_deployed = start_date, 
                      date_retrieved = end_date, 
                      MTexcep = MTlux_exception)
  
  if (is.null(lu) | nrow(lu) == 0) {
    message(paste0('Geo (', i, ') ', geo_trip, ' has no light data - Skipped'))
    GLS.files2[i, c(21:23)] <- list(NA, 'No data', NA)
    next
  }
  
  # Filter to remove noise light during the night
  #d.lux$Light <- lightFilter(d.lux$Light, iter = 2)# 2 iters are enough
  
  #Determine the recording interval
  li<-median(difftime(lu$dtime[2:nrow(lu)],lu$dtime[1:(nrow(lu)-1)],units='mins'))
  
  # Calculate the start and end year
  start_year <- year(lu[1, 1])
  end_year <- year(lu[nrow(lu), 1])
  
  # Check if the GLS trip contains enough data for analysis
  if (start_year == end_year) {
    # Not enough data, indicate for removal
    GLS.files2[i, c(21:23)] <- list(NA, 'Data', 'Remove')
    message(paste0('Remove ', geo_trip, ' without enough data'))
  }
  

  # Select the equinoxes of this year
  Sept_eq <- pull(equinox[equinox$Year == start_year,2])
  Mar_eq <- pull(equinox[equinox$Year == end_year,4])


  ############################################################################################################
  
  # 4. Automated calculation of twilight events ####
  
  # Better use logarithmic values to calculate .trn by twilightCalc
  #lu$lux2<-log(lu$lux)
  #head(lu)
  
  # This function calculates automatically twl for all the period of light data
  twl <- try(
    twilightCalc(lu$dtime, lu$lux, 
                 ask = FALSE, preSelection = TRUE, 
                 LightThreshold = LightThreshold, maxLight = li),
    silent = TRUE
  )
  
  if (inherits(twl, 'try-error')) {
    message(paste('Geo (', i, ') ', geo_trip,' wihtout twilight calculation - Skipped'))
    GLS.files2[i, c(21:23)] <- list(NA, "Data. Problem with twilight calculation", NA)
    next
  }

  #twl<-distinct(twl)
  # 4.1. Filter unrealistic twilight events ####

  #type=1 and type=2 not always correctly assigned, this fixes most of these cases:
  for(j in 1:(nrow(twl)-2)){
    if(sum(twl$type[j:(j+2)])==3){twl$type[j+1]<-2}
    if(sum(twl$type[j:(j+2)])==6){twl$type[j+1]<-1}}
  
  #twl<-distinct(twl)
  
  #Geolight is removing a minute extra unlike bastrack and intiproc when advancing sunset, 
  #I add a minute to compensate
  
  twl$tFirst[twl$type == 2] <- twl$tFirst[twl$type == 2] + 60  
  twl$tSecond[twl$type == 1] <- twl$tSecond[twl$type ==   1] + 60
  
  twl<-twl[,c(1:3)]
  
  # Open a pdf to record the plots that illustrate the edition of twilights
  # pdf(file= paste0(wd,'/Plots/twlEdit/', geo_trip,'_', year_rec, '_', colony,'.pdf'))
  par(mfrow= c(2,1))
  
  # 4.2. Remove and clean unrealistic twilight ####
  twl2<-twilight_cleanup(df = twl,
                         breedingloc_lon = col_lon,
                         breedingloc_lat = col_lat,
                         #days_breeding = c(100:300), #define months when the animal is breeding in the colony
                         months_breeding = c(5,6,7,8), #define months when the animal is breeding in the colony
                         months_extensive_legtucking = NA, 
                         show_plot = T)
  twl2<-distinct(twl2)
  
  twl3<-move_twilights(df = twl2,
                       minutes_different= 15,
                       sun = sun_angle,
                       show_plot = TRUE)
  twl3<-distinct(twl3)
  
  # 4.4. Filter by loess distribution and distanceFilter.

  # We only use the filter of extreme speeds
  twl3<-twl3%>%
    mutate(loess = loessFilter(twl3, k=3, plot = T),
           speed_filter = distanceFilter(twl3, degElevation = sun_angle, distance = 500, 
                                         units = "hour"))
    #filter(loess == TRUE)%>%
    #filter(speed_filter == TRUE)
  
  # dev.off()
  # Adding some columns to read it easily by any algorithm functions
  twl4<-twl3%>%
    mutate(Twilight = tFirst,
           Rise = ifelse(type=='1',TRUE,FALSE),
           confidence = 9,
           Julian = yday(Twilight))
  
  ############################################################################################################
  
  # 5. Calculate raw and smoothed positions ####
  
  # Preliminary sight of the results of the calculations of twilights
  # We smooth in two ways
  # Also, we assign the equinox periods according to the breeding colony latitude
  twl4<-double_smoothing(df=twl4,sun=sun_angle)
  twl4$eq_filter<-assign_equinox_periods(twl4, breedingloc_lat=col_lat)
  
  twl4B <- twl4%>% filter(speed_filter==TRUE)
  # Open a pdf to map raw and smoothed positions
  # pdf(file= paste0(wd,'/Plots/twlEdit/', geo_trip,'_', year_rec, '_', colony,'_rawcoords_smooth.pdf'))
  par(mfrow= c(1,2))


  # Raw
  tripMap(cbind(twl4B$lon_smooth2[twl4B$eq_filter==1],twl4B$lat_smooth2[twl4B$eq_filter==1]), ylim=c(-50,50),xlim=c(-90,50),equinox=F,col='black',legend=T, cex=0.2,add=F) 
  points(cbind(twl4B$lon_smooth2[twl4B$eq_filter==1],twl4B$lat_smooth2[twl4B$eq_filter==1]),cex=0.6,pch=16,col='firebrick')
  points(cbind(col_lon,col_lat),cex=2,pch=16,col='orange')
  
  # Smoothed
  months_breeding=c(5,6,7,8)
  
  tripMap(cbind(twl4B$lon_smooth2[twl4B$eq_filter==1 & !(month(twl4B$tFirst)%in%months_breeding)],twl4B$lat_smooth2[twl4B$eq_filter==1 & !(month(twl4B$tFirst)%in%months_breeding)]), ylim=c(-50,50),xlim=c(-90,50),equinox=F,col='black',legend=T, cex=0.2,add=F) 
  points(cbind(twl4B$lon_smooth2[twl4B$eq_filter==1 & !(month(twl4B$tFirst)%in%months_breeding)],twl4B$lat_smooth2[twl4B$eq_filter==1 & !(month(twl4B$tFirst)%in%months_breeding)]),cex=0.6,pch=16,col='steelblue')
  points(cbind(col_lon,col_lat),cex=2,pch=16,col='orange')
  
  # dev.off()
  # Assign confidence for equinox periods as we used to do 
  twl4 <- twl4 %>%
    mutate(
      confidence = case_when(
        between(Julian, Sept_eq - 20, Sept_eq + 20) ~ 2,
        between(Julian, Sept_eq - 10, Sept_eq + 10) ~ 1,
        between(Julian, Mar_eq - 20, Mar_eq + 20) ~ 2,
        between(Julian, Mar_eq - 10, Mar_eq + 10) ~ 1,
        TRUE ~ confidence
      )
    )

  
 
  ############################################################################################################
  
  # 6. Saving trn files ####
  
  # Select columns to save
  trn<- twl4[, c(1:3,6,7,9,8,18,4,5,16,17)]
  # trn <- distinct(trn)
  
  # Name of the file
  file.trn<-paste0(geo_trip,'_', info$Sp,'_', year_rec,'_',info$Colony, '_automated.trn')

  write.table(trn, 
              paste0(wd, '/Raw Data/GLS/trn/trn_twlcalc/', file.trn), 
              sep = ',', col.names =  T, row.names = F, quote = F)
  
  
  # Saving a trn file in old format to see the results in Locator (BirdTracker)
  trn.bak <- trn %>%
    mutate(
      date_time = format(Twilight, format = '%d/%m/%y %H:%M:%S'),
      type = ifelse(Rise == T, 'Sunrise', 'Sunset')
    ) %>%
    arrange(date_time) %>%
    dplyr::select(date_time, type, confidence)
 
  # Change the format of the date column to match with 'Locator' software
  
  file.trn.bak<-paste0(geo_trip,'_', info$Sp,'_', year_rec,'_',info$Colony, '.trn')
  
  #Save trn file
  write.table(trn.bak, 
              paste0(wd, '/Raw Data/GLS/trn/trn_backup/', file.trn.bak), 
              sep = ',', col.names =  F, row.names = F, quote = F)
  
  
 
  # Add the name of the file to the dataset
  cat('Trn from geo (', i,') ',geo_trip, ' is saved \n')
  GLS.files2[i, c(21:23)] <- list(file.trn, 'Data', 'No')
 
  trn.final1<-cbind(geo_trip,trn)
  trn.final1<-as.data.frame(trn.final1)
  trn.list1[i] <- list(trn.final1)
  trn.final2<-cbind(geo_trip,trn.bak)
  trn.final2<-as.data.frame(trn.final2)
  trn.list2[i] <- list(trn.final2)
}
time<-toc()# it spent 2 hours 25 minutes to run almost 1300 trips
GLS.files2
#GLS.files2_bak <- GLS.files2
#### 3. Save files and trn data ####

# Filtering trips we need to do again
GLS.files2<-GLS.files2%>%
  filter(Do_again != 'Remove' | is.na(Do_again))

write_csv(GLS.files2,paste0(wd, '/Raw Data/GLS.files2.csv'))

# GLS.files3<-GLS.files2%>%
#   filter(Do_again == 'Yes' | Do_again == 'Remove')
#write.csv(GLS.files3,paste0(wd, '/Raw Data/GLS.files_doagain.csv')) # New .csv dataset contains the files we need to do again with new end dates


trn.all1 <- ldply(trn.list1, data.frame)
names(trn.all1)
head(trn.all1)
colnames(trn.all1)[1]<-'geo_trip'
trn.all1$geo_trip<-as.character(trn.all1$geo_trip)

trn.all2 <- ldply(trn.list2, data.frame)
names(trn.all2)
head(trn.all2)
colnames(trn.all2)[1]<-'geo_trip'
trn.all2$geo_trip<-as.character(trn.all2$geo_trip)

write_csv(trn.all1, paste0(wd, '/Raw Data/trn_all.csv'))# save as csv

