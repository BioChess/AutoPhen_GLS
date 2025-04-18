#### SET UP GEOLOCATOR ANALYSES ####

#### Updated: 23-10-2023

#### Created by Diego Vicente-Sastre
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals, 
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)


################################################################################

# Function to read any light data file from different models producers of 
# light-level geolocators
# @producer : brand or producer of the model
# @geo_trip : ID of the geolocator or trip
# @file_name : name of the light file to open
# @start_date : starting date the data to analyse
# @end_Date : ending date the data to analyse
# @MTexcep : IDs of those loggers from Migrate Technology, that cannot be opened
# using the readMTlux function

read_lightGLS <- function(producer, ID, file, 
                          start = NULL, end = NULL, MTexcep = NULL) {
  require(lubridate)
  if (file.exists(file)) {
    if (Producer %in% c('BAS', 'Biotrack', 'Lotek')) {
      lu <- read.table(file, sep = ",", skip = 1, header = FALSE, fill = TRUE)
      lu$Date <- parse_date_time(lu[, 2], 
                                 orders = c("dmy HMS", "dmY HMS",
                                            "ymd HMS", "mdy HMS"), 
                                 tz = "GMT")
      lu$Light <- lu[, 4]
      lu <- lu[, c(5, 6)]
    } else if (Producer == 'Migrate Technology') {
      if(ID %in% MTexcep) {
        cat(paste0('File from geo ', geoID, ' cannot be opened by readMTlux function'))
        lu<-read.csv(file, skip = 19, sep = '\t', 
                     col.names = c('Date', 'Light'), 
                     colClasses = c('character','numeric'))
        lu$Date <- as.POSIXct(strptime(lu$Date, "%d/%m/%Y %H:%M:%S", tz = "GMT"))
        lu<-lu%>%dplyr::filter(!is.na(Date))
      } else{lu<-readMTlux(file)}
    }
    
    lu <- lu %>%
      mutate(dtime = as.POSIXct(Date, tz = 'GMT'),
             lux = as.numeric(gsub('\\,', '.', Light)),
             date = as.Date(Date),
             time = strftime(Date, format = "%H:%M:%S"))
    if(!is.null(start) & !is.null(end)){
      lu <- lu %>%
        filter(date >= start & date <= end)
    }

    
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

calib_period<- function (data, start, end, threshold, lon, lat) {
  
  PD <- data %>% filter(Date > start & Date < end)
  
  if (nrow(PD) == 0 || round(as.numeric(difftime(max(PD$Date), min(PD$Date), 
                                                 units = 'days')))*2 < 3) {
    cat('No calibration for step', step, '\n')
    return(c('NA', 'NA', 'NA', 'NA'))
  } else {
    lightImage(tagdata = PD, offset = 12, zlim = c(0, 20))
    tsimageDeploymentLines(PD$Date, lon = lon, lat = lat, offset = 12, lwd = 3, 
                           col = adjustcolor('orange', alpha.f = 0.5))
    tm_step <- as.POSIXct(c(start, end), tz = 'GMT')
    abline(v = tm_step, lwd = 2, lty = 2, col = 'orange')
    twl <- BAStag::preprocessLight(PD, threshold = threshold, offset = 12, 
                                   lmax = 12, map = TRUE, plotMap = plotMapp)
    twl <- subset(twl, Deleted == 'FALSE')
    return(twl)
  }
}


################################################################################

# Function to check and replace outliers with NA values
# After all loggers are calibrated and twilight errors are calculated,
# we check possible outliers and uncoherent values
# @data : summary of twilight errors and zenith angles of each logger
# @var_name : vector of variables' names to calculate and replace outliers

calib_out <- function(data, vars) {
  for (v in vars) {
    # v = vars[2]
    if(all(is.na(data[[v]]))) {
      data
    } else{
         b <- boxplot(data[[v]])
    out <- b$out
    data <- data %>%
      mutate({{v}} := ifelse(data[[v]] %in% out, NA, data[[v]])
      )
  }
    }
  return(data)
}

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

select_sunAngle <- function(producer, colony, year, calib.sum, calibs, param = T) {
  
  
  calib <- calibs %>%
    filter(geo_alfanum == geoID)
  if (nrow(calib) == 0) {
    cat(paste0('The logger has no calibrations \n'))
    calib <- calib.sum %>%
      filter(Producer == producer, Colony == colony, Year == year)
  } else {
    cat(paste0('The logger has calibrations \n'))
  }
  
  if (nrow(calib) > 0) {
    sun_angle <- 90 - calib$mean_zenith
    calib.param <- list(zenith_mean = calib$mean_zenith, 
                        zenith0_mean = calib$mean_zenith0, 
                        alpha =  c(calib$mean_alpha.mean, calib$mean_alpha.sd))
  } else {
    message('No calibration data found. Sun_angle set to NA \n')
    sun_angle <- NA
  }
  
  if(param == F){return(sun_angle)}
  if(param == T){return(calib.param)}
  
}

################################################################################

calc_tol <- function(twl, n = 0.001, Sep_eq, Mar_eq, d, k1 = 0.95, k2= 0.5, plot = T){
  
  require(SGAT)
  require(patchwork)
  
  tsep <- which(as.Date(twl$Twilight) >= (Sep_eq - d) & 
                  as.Date(twl$Twilight) <= (Sep_eq+ d))
  tmar <- which(as.Date(twl$Twilight) >= (Mar_eq - d) & 
                  as.Date(twl$Twilight) <= (Mar_eq+ d))
  t <- c(tsep, tmar)
  
  # # Define 10.000 posible tols
  tol_seq <-seq(from = n, to=0.2, by = n)
  
  # Create an empty matrix to store the results
  m <- matrix(0, nrow = length(tol_seq), ncol = 6)
  colnames(m) <- c('tol', 'dif.lat', 'dif.lat.sep', 'dif.lat.mar', 
                               'dif.lon.sep', 'dif.lon.mar')
  
  for (i in seq_along(tol_seq)) {
    # j = 1
    tol <- tol_seq[i]
    # print(tol)
    
    path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = tol)
    x0 <- path$x
    
    xsep <- x0[tsep, ]
    xmar <- x0[tmar, ]
    # mean(xmar)
    x <- x0[t, ]
    
    max.lat.sep <- max(xsep[, 2])
    min.lat.sep <- min(xsep[, 2])
    dif.lat.sep <- max.lat.sep - min.lat.sep
    
    max.lon.sep <- max(xsep[, 1])
    min.lon.sep <- min(xsep[, 1])
    dif.lon.sep <- max.lon.sep - min.lon.sep
    
    if(is.matrix(xmar) == T){
      max.lat.mar <- max(xmar[, 2])
      min.lat.mar <- min(xmar[, 2]) 
      dif.lat.mar <- max.lat.mar - min.lat.mar
      
      max.lon.mar <- max(xmar[, 1])
      min.lon.mar <- min(xmar[, 1])
      dif.lon.mar <- max.lon.mar - min.lon.mar
      
    } else if (is.matrix(xmar) == F){
      dif.lat.mar <- Inf
    }
    
    max.lat <- max(x[, 2])
    min.lat <- min(x[, 2])
    dif.lat <- max.lat - min.lat
    
    m[i, ] <- c(tol, dif.lat, dif.lat.sep, dif.lat.mar,
                            dif.lon.sep, dif.lon.mar)
  }
  
  # Convert the result_matrix into a data frame
  tol_df <- as.data.frame(m)
  
  tol_df <- tol_df %>%
    mutate(across(c(dif.lat, dif.lat.sep, dif.lat.mar), round, digits = 2))%>%
    mutate(dec_sep = lag(dif.lat.sep) - dif.lat.sep,
           dec_mar = lag(dif.lat.mar) - dif.lat.mar) 
  
  tols_df <- tol_df %>%
    filter(dec_sep!= 0) %>%
    mutate(cumsum_sep = cumsum(replace_na(dec_sep, 0)))
  tolm_df <- tol_df %>%
    filter(dec_mar!= 0) %>%
    mutate(cumsum_mar = cumsum(replace_na(dec_mar, 0)))
  
  
  tol_s <- tols_df$tol[which(tols_df$cumsum_sep >= max(tols_df$cumsum_sep)*k1)[1]]
  
  if(plot == T){
    tols_plot <- ggplot(tols_df)+
      geom_point(aes(x = tol, y = cumsum_sep), color = 'grey25', alpha = 0.25)+
      geom_smooth(aes(x = tol, y = cumsum_sep), 
                  method = "lm", formula = y ~ log(x), se = FALSE)+
      geom_vline(xintercept = tol_s, linetype = "dashed", color = "red", linewidth = 1.2)+ 
      theme_bw()+
      theme(plot.title = element_text(size = 20),
            axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_text(size = 16),
            axis.text.y = element_text(size = 14))+
      labs(title = 'Latitude difference along fall equinox')+
      ylab('Latitude difference') + xlab('tol value')
  }

  
  if(length(tmar) >= 15){
    
    tol_m <- tolm_df$tol[which(tolm_df$cumsum_mar >= max(tolm_df$cumsum_mar)*k2)[1]]
    if(plot == T){
      tolm_plot <- ggplot(tolm_df)+
        geom_point(aes(x = tol, y = cumsum_mar), color = 'grey25', alpha = 0.25)+
        geom_smooth(aes(x = tol, y = cumsum_mar),
                    method = "lm", formula = y ~ log(x), se = FALSE)+
        geom_vline(xintercept = tol_m, linetype = "dashed", color = "red", linewidth = 1.2)+ 
        theme_bw()+
        theme(plot.title = element_text(size = 20),
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))+
        labs(title = 'Latitude difference along spring equinox')+
        ylab('Latitude difference') + xlab('tol value')
    }
  } else {
    tol_m <- NA
    if(plot == T){
      tolm_plot <- ggplot(tol_df)+
        labs(title = 'Latitude difference along spring equinox',
             xlab ='tol value',
             ylab = ' Latitude difference')+
        xlim(c(0.00, 0.20))+ ylim(c(0, 100))+
        theme_bw()+
        theme(plot.title = element_text(size = 20),
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
    }

  }
  
  cat(paste0('Tol values calculated'),'\n')
  
  
  if(plot == T){
    p <-plot_grid(tols_plot / tolm_plot)
  } else{
    p <- NULL
  }
  
  out <- list(tol_sep = tol_s, tol_mar = tol_m, plot = p)
  return(out)
  
}

################################################################################

plot_SGAT <- function (x0, time, d, Sep_eq, Mar_eq, world, col_lat, col_lon){
  
  require(rnaturalearth)
  require(ggplot2)
  
  
  xsep <- x0[which(as.Date(time) >= Sep_eq - d & 
               as.Date(time) <= Sep_eq + d),]
  xmar <- x0[which(as.Date(time) >= Mar_eq - d & 
                     as.Date(time) <= Mar_eq + d),]
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  p <- ggplot()+
    geom_sf(data = world, fill=grey(0.9), color=grey(0.6), lwd = 0.2)+
    coord_sf(xlim = c(min(x0[,1]), max(x0[,1])),
             ylim = c(min(x0[,2]), max(x0[,2])), expand = T) +
    geom_path(data = data.frame(x0), aes(x = lon, y = lat), color = "firebrick") +
    geom_point(data = data.frame(x0), aes(x = lon, y = lat), shape = 19, color = "firebrick") +
    geom_point(data = data.frame(xsep), aes(x = lon, y = lat), shape = 19, color = "cornflowerblue") +
    geom_point(data = data.frame(xmar), aes(x = lon, y = lat), shape = 19, color = "darkgreen") +
    geom_point(aes(x = col_lon, y = col_lat), shape = 16, size = 2.5, color = "black") +
    xlab('Longitude') + ylab('Latitude')+
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14))
  
  return(p)
}
