## FIX RESUSCITATED TRN FILES ##

# Script developed by Diego Vicente-Sastre
# Based on previous script developed by SEATRACK project

# Updated: 05-03-2023


############################################################################################################

#### 0.1 Working directory ####
Sys.setenv(TZ = 'GMT') 
### !!! IMPORTANT: IF NOT SET LIKE THIS, MAKE PROBLEMS TO CREATE DATE_TIME FROM 
# PASTING DATE plus TIME

place <- 'PCUB'  # pc in the UB
#place <- "home"  # laptop at home


if (place == 'PCUB') WD <- 'D:/Dropbox'
if (place == 'home') WD <- 'F:/Dropbox'  # pc in the UB
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
               'tictoc', 'cowplot', 'caTools')
detach("package:MASS", unload = TRUE)
detach("package:mice", unload = TRUE)
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
GLS.files <- read_csv('Raw Data/GLS.files2.csv', col_names = TRUE)
GLS.info <- loadGLSInfo() 

colonies <- loadCol.Coord() # colony's coordinates

calib.info <- loadCalib.Info() # calibration periods
calibrations <- loadCalib.Angles() # calibration parameters calculated
calib.summary <- loadCalib.Summary() # calibration summary per year and colony

#Equinox vary a little between years, we selecte the precise dates for each GLS
equinox<-read_excel('Raw Data/Equinox_by_year.xlsx')
head(equinox)
colony_days<-read_csv('Raw Data/colony_days_v3.csv')
# Define threshold and offset
offset = 12 # adjusts the y-axis to put night (dark shades) in the middle 
LightThreshold <- 1.5

############################################################################################################

# List of file paths for resuscitated files
file_paths <- c(
  # file.path(WD, 'Alumnos/Sara Rubio/Revision_resucitats_Sara_DV.xlsx'),
  file.path(WD, 'Alumnos/Revision_resucitats_The_DV.xlsx'),
  file.path(WD, 'Alumnos/Mireia_Marrugat/Datasets/Geos/resucitats/Revisio_resucitats_Mireia.xlsx'),
  file.path(WD, 'Alumnos/Revision_resucitats_Jul_DV.xlsx')
)

# Load and filter resuscitated files
res_alum <- loadAndFilterResuscitated(file_paths)

bad_trips <- read_csv('Raw Data/bad_trips.csv')

bad_GLS.files <- GLS.files %>%
  left_join(bad_trips, .)

GLS.files2 <- GLS.files %>%
  filter(!(geo_trip %in% bad_trips[bad_trips$Problem == 'Removed',]$geo_trip))

light_path <- paste0(wd,'/Raw Data/GLS/lux/')
LightThreshold<-1.5
MTlux_exception <- c('N643003001', 'Q666002001', 'X250001001')

trn_path <- paste0(wd, '/Raw Data/GLS/trn/trn_twlcalc/')
trn.list <-list()
for (i in 1:length(bad_GLS.files)) {
  #i= 96
  print(i)

  # 2. Extract the GLS information ####
  info<-bad_GLS.files[i,]
  
  lightfile <- paste0(light_path,if_else(!is.na(info$File.adjlux), info$File.adjlux, info$File.lux))
  trn_file <- paste0(trn_path, info$File.trn)
  trn<-read_csv(trn_file, col_names = T, show_col_types = FALSE)
  geo_trip<-info$geo_trip
  geoID <- info$geo_alfanum
  year_rec <- info$Year_rec #which year we recovered the GLS
  col.name <- info$Colony #which colony is from
  if(col.name == 'M.Clara') col.name <- 'MClara'
  if(col.name == 'Cala Morell') col.name <- 'CalaMorell'
  if(col.name == 'Curral Velho') col.name <- 'CVelho'
  if(col.name == 'Ilheu de Cima') col.name <- 'ICima'
  colony<- colonies%>%
    filter(Colony == col.name)%>%
    dplyr::select(Longitude = Lon, Latitude = Lat)
  
  Producer<- info$Producer #company
  
  start_date<-as.POSIXct(info$start_date, tz = 'GMT') #date start the trip
  end_date <-as.POSIXct(info$end_date, tz = 'GMT') #date end the trip
  
  cdays<-colony_days%>%
    dplyr::filter(geo_trip == .env$geo_trip)%>%
    filter(light == T | wetdry == T)
  
  problem <- info$Problem
  
  if(problem == 'Removed'){
    cat('Geo (', i,') ', geo_trip, 'for year', year_rec, ' from ', col.name,  ' has been removed\n')
    next
  }
  
  cat('Processing geo (', i,') ', geo_trip, 'for year', year_rec, ' from ', col.name,  '\n')
  #### 2.2. Read light data ####  
  
  lu <- read_lightGLS(Producer, geo_trip, 
                      file_name = lightfile, 
                      date_deployed = start_date, 
                      date_retrieved = end_date, 
                      MTexcep = MTlux_exception)
  
  #Filter dates between deployment and recovery
  lu2 <- lu %>%
    mutate(date = as.Date(Date)) %>%
    filter(date >= start_date & date <= end_date)
    
    
  if(problem == 'Cut' | problem == 'Incomplete'){  
    # Days without light
    nolu <- lu2 %>%
      group_by(date) %>%
      mutate(Light =  ifelse(Light > 64, 64, Light))%>% # adapt to Biotrack format
      dplyr::summarise(sum.Light = sum(Light)) %>%
      mutate(percent = round(sum.Light/max(sum.Light) * 100, digits = 2))%>%
      bind_cols(geo_trip = geo_trip, Producer = Producer)%>%
      relocate(geo_trip, .before = 'date') %>%
      relocate(Producer, .after = 'geo_trip')%>%
      mutate(diff_percent = abs(percent -lag(percent)),
             diff_percent2 = abs(runmin(percent -lag(percent), k = 3)))
    
    sd_diff <- sd(nolu$diff_percent2, na.rm = T)
    mean_diff <- mean(nolu$diff_percent2, na.rm = T)
    
    date_to_cut <- nolu %>%
      filter(percent == 100 & diff_percent2 > mean_diff+sd_diff)%>%
      pull(date)
    
    # if(length(date_to_cut) == 0){
    #   date_to_cut <- nolu %>%
    #     filter(diff_percent2 > mean_diff+sd_diff)%>%
    #     pull(date)
    # }
    
    date_to_cut <- min(date_to_cut)
    
    nolu <- nolu %>%
      mutate(color = ifelse(date >= date_to_cut, 'firebrick', 'steelblue'))
    ggplot(nolu, aes(x=date, y= percent, color = color))+
      geom_point()+
      scale_color_identity(guide = 'legend',
                           name = 'Dates to be removed',
                           breaks = c('firebrick', 'steelblue'),
                           labels = c('Dates removed', 'Dates preserved')
      )+
      geom_vline(xintercept = date_to_cut, linetype = 'dashed')+
      ylim(0,100)+
      labs(title = geo_trip,
           x = 'Date',
           y = 'Daily light intensity (%)')+
      theme_minimal()
    
    ggplot(nolu, aes(x=date, y= diff_percent2, color = color))+
      geom_point()+
      scale_color_identity(guide = 'legend',
                           name = 'Dates to be removed',
                           breaks = c('firebrick', 'steelblue'),
                           labels = c('Dates removed', 'Dates preserved')
      )+
      geom_vline(xintercept = date_to_cut, linetype = 'dashed')+
      ylim(-100,100)+
      labs(title = geo_trip,
           x = 'Date',
           y = 'Daily light intensity (%)')+
      theme_minimal()
    
    
    lu3 <- lu2 %>%
      filter(date <= date_to_cut)
    
  }
  

  if(problem == 'Moved'){
    next
    sun_angle <- select_sun_angle(producer = Producer, col.name, year = year_rec,
                                  calib.summary, calibrations, param = F)
    
    li<-median(difftime(lu$dtime[2:nrow(lu)],lu$dtime[1:(nrow(lu)-1)],units='mins'))
    twl <- try(
      twilightCalc(lu2$dtime, lu2$lux,
                   ask = FALSE, preSelection = TRUE,
                   LightThreshold = LightThreshold, maxLight = li),
      silent = TRUE
    )

    
    twl_clean<-twilight_cleanup(df = twl,
                                breedingloc_lon = colony$Longitude,
                                breedingloc_lat = colony$Latitude,
                                # days_breeding = c(100:300), #define months when the animal is breeding in the colony
                                months_breeding = c(5,6,7,8), #define months when the animal is breeding in the colony
                                months_extensive_legtucking = NA,
                                show_plot = T)
    
    twl_removed<-move_twilights(df = twl_clean,
                                minutes_different= 15,
                                sun = sun_angle,
                                show_plot = TRUE)
    
    
    trn <- twl_removed%>%
      mutate(Twilight = tFirst,
             Rise = ifelse(type=='1',TRUE,FALSE),
             confidence = 9,
             Julian = yday(Twilight))
    
    trn <- trn %>%
      filter(!(as.Date(Twilight) %in% cdays$date))
    # Find time differences in twilight patterns
    # trn <- trn[-c(1:2, (nrow(trn)-10):nrow(trn)),]
    # Sunrises
    sunr <- trn %>%
      filter(Rise == T)%>%
      mutate(dec_time = round(hour(Twilight)+minute(Twilight)/60+second(Twilight)/3600, digits = 2))%>%
      mutate(diff_time = round(c(NA, diff(dec_time)), digits = 2)) %>% 
      mutate(diff_time = ifelse(abs(diff_time) > 20, NA, diff_time))
    
    sunr_row <- which.max(abs(sunr$diff_time))
    sunr_change <- as.Date(sunr$Twilight[sunr_row])
    sunr_diff <- sunr$diff_time[sunr_row]
    
    # Sunsets
    suns <- trn %>%
      filter(Rise == F)%>%
      mutate(dec_time = round(hour(Twilight)+minute(Twilight)/60+second(Twilight)/3600, digits = 2))%>%
      mutate(diff_time = round(c(NA, diff(dec_time)), digits = 2))%>% 
      mutate(diff_time = ifelse(abs(diff_time) > 20, NA, diff_time))
    
    suns_row <- which.max(abs(suns$diff_time))
    suns_change <-  as.Date(suns$Twilight[suns_row])
    suns_diff <- suns$diff_time[suns_row]
    
    row_change <- c(sunr_row, suns_row)[which.max(abs(c(sunr_diff, suns_diff)))]
    twl_diff <- c(sunr_diff, suns_diff)[which.max(abs(c(sunr_diff, suns_diff)))]
    
    if(geo_trip == '1643001001'){
      suns_diff <- -16.75
      sunr_diff <- -16.75
      sunr_row <- 1
      suns_row <- 1
    } 
    if(geo_trip == '19426003001'){
      suns_diff <- suns_diff-0.5
      sunr_diff <- suns_diff
    }
    if(geo_trip == '19426003002'){
      suns_diff <- suns_diff-3.30
      sunr_diff <- suns_diff
    }
    if(geo_trip == '19922003001'){
      suns_diff <- 16.92-20.18 # there is an outlier in the middle. Removed, plot it is not good
    }
    if(geo_trip == '19959003001'){
      sunr_row <- suns_row
      suns_diff <- suns_diff-4
      sunr_diff <- sunr$diff_time[sunr_row]+2
    }
    if(geo_trip == '21166002001'){
      suns_diff <- -3
      sunr_diff <- -3
      sunr_row <- 1
      suns_row <- 1
    }
    if(geo_trip == '21166002001'){
      suns_diff <- -3
      sunr_diff <- -3
    
    }
    
    if(geo_trip == '23297004001'){
      suns_diff <- 0.3
      sunr_diff <- 0.3
      
    }
    
    if(geo_trip == '23348004001'){
      suns_diff <-  sunr_diff+0.3
    }
    
    if(geo_trip == '23355004001'){
      sunr_diff <- 5
      suns_diff <- 5
    }
    
    # Correct sunrises
    
    sunr2 <- sunr %>%
      mutate(Modify = ifelse(row_number() < sunr_row, F, T))%>%
      mutate(Twilight2 = if_else(Modify == F, Twilight, Twilight-(sunr_diff*3600)))%>%
      mutate(dec_time2 = hour(Twilight2)+minute(Twilight2)/60+second(Twilight2)/3600)
    
    # Correct sunsets
  
    suns2 <- suns %>%
      mutate(Modify = ifelse(row_number() < suns_row, F, T))%>%
      mutate(Twilight2 = if_else(Modify == F, Twilight, Twilight-(suns_diff*3600)))%>%
      mutate(dec_time2 = hour(Twilight2)+minute(Twilight2)/60+second(Twilight2)/3600)
    
    # suns_diff2 <- suns$diff_time2[which.max(abs(suns$diff_time2))]
    # plot(as.Date(sunr$Twilight), sunr$diff_time)
    # plot(as.Date(suns$Twilight2), suns$dec_time2)
    
    trn2 <- rbind(sunr2, suns2)
    trn2 <- trn2 %>% arrange(Twilight2)
    
    ggplot()+
      geom_point(data = trn2, aes(x = as.Date(Twilight), y = dec_time, color = type))+
      geom_vline(xintercept = suns_change)+
      geom_vline(xintercept = sunr_change, linetype = 'dashed')
    ggplot()+
      geom_point(data = trn2, aes(x = as.Date(Twilight2), y = dec_time2, color = type))+
      geom_vline(xintercept = suns_change)+
      geom_vline(xintercept = sunr_change, linetype = 'dashed')
  }
  
  # path <- thresholdPath(trn2$Twilight2, trn2$Rise, zenith = sun_angle-90, tol=0.18)
  # x0 <- path$x
  # raw_map <- ggplot()+
  #   geom_sf(data = world, fill=grey(0.9), color=grey(0.6), lwd = 0.2)+
  #   coord_sf(xlim = c(min(x0[,1]), max(x0[,1])),
  #            ylim = c(min(x0[,2]), max(x0[,2])), expand = T) +
  #   geom_path(data = data.frame(x0), aes(x = lon, y = lat), color = "firebrick") +
  #   geom_point(data = data.frame(x0), aes(x = lon, y = lat), shape = 19, color = "firebrick") +
  #   labs(title = geo_trip)+
  #   xlab('Longitude') + ylab('Latitude')+
  #   theme_minimal() +
  #   theme(panel.grid = element_blank(),
  #         panel.border = element_rect(color = "black", fill = NA),
  #         plot.title = element_text(size = 20),
  #         axis.title.x = element_text(size = 16),
  #         axis.text.x = element_text(size = 14),
  #         axis.title.y = element_text(size = 16),
  #         axis.text.y = element_text(size = 14))
  # raw_map

  if(Producer == 'Biotrack') end  <- '.lig'
  if(Producer == 'Migrate Technology') end <- '.lux'
  file.ligth<-paste0(geo,'_', info$Sp,'_', year_rec,'_',colony, '_fixed', end)
  write.table(trn2,
              paste0(wd, '/Raw Data/GLS/trn/trn_fixed/', file.trn),
              sep = ',', col.names =  T, row.names = F, quote = F)

}

trn.all <- plyr::ldply(trn.list, data.frame)
write_csv(trn.all, paste0(wd, '/Raw Data/cut_dates_resuscitated.csv'))# save as csv


