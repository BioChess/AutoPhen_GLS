#### SET UP GEOLOCATOR ANALYSES ####

#### Updated: 2024 (corrected version for course)

#### Created by Diego Vicente-Sastre
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals,
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)


################################################################################

# Function to read any light data file from different models/producers of
# light-level geolocators.
#
# @producer   : brand or producer of the model
# @ID         : ID of the geolocator
# @file       : path to the light file
# @start      : (optional) starting date to subset data
# @end        : (optional) ending date to subset data
# @MTexcep    : (optional) IDs of Migrate Technology loggers that cannot be
#               opened with readMTlux

read_lightGLS <- function(producer, ID, file,
                          start = NULL, end = NULL, MTexcep = NULL) {
  require(lubridate)

  if (!file.exists(file)) {
    message(paste0('File from geo ', ID, ' does not exist'))
    return(NULL)
  }

  # FIX: use argument 'producer' (not global 'Producer')
  if (producer %in% c('BAS', 'Biotrack', 'Lotek')) {
    lu <- read.table(file, sep = ",", skip = 1, header = FALSE, fill = TRUE)
    lu$Date  <- parse_date_time(lu[, 2],
                                orders = c("dmy HMS", "dmY HMS",
                                           "ymd HMS", "mdy HMS"),
                                tz = "GMT")
    lu$Light <- lu[, 4]
    lu <- lu[, c(5, 6)]

  } else if (producer == 'Migrate Technology') {
    # FIX: use argument 'ID' (not global 'geoID') in message and condition
    if (!is.null(MTexcep) && ID %in% MTexcep) {
      cat(paste0('File from geo ', ID, ' cannot be opened by readMTlux function\n'))
      lu <- read.csv(file, skip = 19, sep = '\t',
                     col.names   = c('Date', 'Light'),
                     colClasses  = c('character', 'numeric'))
      lu$Date <- as.POSIXct(strptime(lu$Date, "%d/%m/%Y %H:%M:%S", tz = "GMT"))
      lu <- lu %>% dplyr::filter(!is.na(Date))
    } else {
      lu <- readMTlux(file)
    }

  } else {
    message('Unknown producer: ', producer, '. Returning NULL.')
    return(NULL)
  }

  lu <- lu %>%
    mutate(dtime = as.POSIXct(Date, tz = 'GMT'),
           lux   = as.numeric(gsub('\\,', '.', Light)),
           date  = as.Date(Date),
           time  = strftime(Date, format = "%H:%M:%S"))

  if (!is.null(start) && !is.null(end)) {
    lu <- lu %>% filter(date >= start & date <= end)
  }

  return(lu)
}


################################################################################

# Helper function: plot a map inside BAStag::preprocessLight.
# Called automatically by calib_period() — do not call directly.
# Relies on 'map', 'xlim' and 'ylim' defined in the calling environment.
#
# @xlim : longitude limits of the map
# @ylim : latitude limits of the map

plotMapp <- function(xlim, ylim) {
  plot.new()
  sp::plot(map, xlim = xlim, ylim = ylim)
  plot.window(xlim, ylim)
}


################################################################################

# Function to perform the interactive calibration for one calibration period.
# Displays the light image and calls preprocessLight for manual twilight review.
# Returns the cleaned twilight data frame, or a vector of NA_real_ if the
# period is too short or contains no data.
#
# @data      : light data frame from read_lightGLS
# @start     : starting POSIXct of the calibration period
# @end       : ending POSIXct of the calibration period
# @threshold : light threshold for preprocessLight
# @lon       : longitude of the calibration site
# @lat       : latitude of the calibration site

calib_period <- function(data, start, end, threshold, lon, lat) {

  PD <- data %>% filter(Date > start & Date < end)

  # FIX: return NA_real_ (not string 'NA') so downstream numeric ops work
  if (nrow(PD) == 0 ||
      round(as.numeric(difftime(max(PD$Date), min(PD$Date),
                                units = 'days'))) * 2 < 3) {
    message('No calibration data or period too short')
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }

  lightImage(tagdata = PD, offset = 12, zlim = c(0, 20))
  tsimageDeploymentLines(PD$Date, lon = lon, lat = lat, offset = 12,
                         lwd = 3, col = adjustcolor('orange', alpha.f = 0.5))
  tm_step <- as.POSIXct(c(start, end), tz = 'GMT')
  abline(v = tm_step, lwd = 2, lty = 2, col = 'orange')

  twl <- BAStag::preprocessLight(PD, threshold = threshold, offset = 12,
                                  lmax = 12, map = TRUE, plotMap = plotMapp)
  twl <- subset(twl, Deleted == FALSE)
  return(twl)
}


################################################################################

# Function to check and replace outliers with NA.
# Applied after all loggers are calibrated to catch incoherent zenith/alpha values.
#
# @data : data frame with calibration parameters (one row per GLS)
# @vars : character vector of column names to check for outliers

calib_out <- function(data, vars) {
  for (v in vars) {
    if (all(is.na(data[[v]]))) next  # FIX: use next instead of bare 'data'
    b   <- boxplot(data[[v]], plot = FALSE)
    out <- b$out
    data <- data %>%
      mutate({{ v }} := ifelse(data[[v]] %in% out, NA, data[[v]]))
  }
  return(data)
}


################################################################################

# Function to select the zenith angle (and optionally alpha parameters) for a
# given GLS from the individual calibration results or, if unavailable, from
# the colony/year/producer summary.
#
# @producer  : producer string (e.g. "Migrate Technology", "Biotrack")
# @geoID     : ID of the geolocator (used to look up individual calibration)
# @colony    : colony name
# @year      : recovery year
# @calib.sum : data frame with summary calibration parameters per colony/year
# @calibs    : data frame with individual calibration parameters per GLS
# @param     : if FALSE return sun_angle scalar; if TRUE return full calib list

select_sunAngle <- function(producer, geoID, colony, year,
                            calib.sum, calibs, param = TRUE) {

  # FIX: geoID is now an explicit argument, not a global variable
  calib <- calibs %>% filter(geo_alfanum == geoID)

  if (nrow(calib) == 0) {
    cat('The logger has no individual calibration — using colony summary\n')
    calib <- calib.sum %>%
      filter(Producer == producer, Colony == colony, Year == year)
  } else {
    cat('The logger has individual calibration\n')
  }

  if (nrow(calib) > 0) {
    sun_angle  <- 90 - calib$mean_zenith
    calib.param <- list(
      zenith_mean  = calib$mean_zenith,
      zenith0_mean = calib$mean_zenith0,
      alpha        = c(calib$mean_alpha.mean, calib$mean_alpha.sd)
    )
  } else {
    message('No calibration data found for this GLS. sun_angle set to NA')
    sun_angle   <- NA
    calib.param <- list(zenith_mean = NA, zenith0_mean = NA, alpha = c(NA, NA))
  }

  if (!param) return(sun_angle)
  if (param)  return(calib.param)
}


################################################################################

# Function to calculate the optimal tolerance (tol) value for thresholdPath,
# by evaluating the latitude spread around equinox periods across a range of
# tol values.
#
# @twl    : twilight data frame with Twilight and Rise columns
# @zenith : zenith angle for thresholdPath
# @n      : step size for tol sequence (default 0.001)
# @Sep_eq : Date of September equinox
# @Mar_eq : Date of March equinox
# @d      : number of days around equinox to evaluate
# @k1     : fraction of cumulative latitude change to select fall tol (default 0.95)
# @k2     : fraction of cumulative latitude change to select spring tol (default 0.5)
# @plot   : whether to display diagnostic plots

calc_tol <- function(twl, zenith, n = 0.001, Sep_eq, Mar_eq, d,
                     k1 = 0.95, k2 = 0.5, plot = TRUE) {

  require(SGAT)
  require(patchwork)

  tsep <- which(as.Date(twl$Twilight) >= (Sep_eq - d) &
                  as.Date(twl$Twilight) <= (Sep_eq + d))
  tmar <- which(as.Date(twl$Twilight) >= (Mar_eq - d) &
                  as.Date(twl$Twilight) <= (Mar_eq + d))
  t <- c(tsep, tmar)

  tol_seq <- seq(from = n, to = 0.2, by = n)

  m <- matrix(0, nrow = length(tol_seq), ncol = 6)
  colnames(m) <- c('tol', 'dif.lat', 'dif.lat.sep', 'dif.lat.mar',
                   'dif.lon.sep', 'dif.lon.mar')

  for (i in seq_along(tol_seq)) {
    tol  <- tol_seq[i]
    # FIX: use argument 'zenith' instead of global variable
    path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = tol)
    x0   <- path$x

    xsep <- x0[tsep, ]
    xmar <- x0[tmar, ]
    x    <- x0[t, ]

    dif.lat.sep <- max(xsep[, 2]) - min(xsep[, 2])
    dif.lon.sep <- max(xsep[, 1]) - min(xsep[, 1])

    if (is.matrix(xmar)) {
      dif.lat.mar <- max(xmar[, 2]) - min(xmar[, 2])
      dif.lon.mar <- max(xmar[, 1]) - min(xmar[, 1])
    } else {
      dif.lat.mar <- Inf
      dif.lon.mar <- Inf
    }

    dif.lat <- max(x[, 2]) - min(x[, 2])

    m[i, ] <- c(tol, dif.lat, dif.lat.sep, dif.lat.mar,
                dif.lon.sep, dif.lon.mar)
  }

  tol_df <- as.data.frame(m) %>%
    mutate(across(c(dif.lat, dif.lat.sep, dif.lat.mar), round, digits = 2),
           dec_sep = lag(dif.lat.sep) - dif.lat.sep,
           dec_mar = lag(dif.lat.mar) - dif.lat.mar)

  tols_df <- tol_df %>%
    filter(dec_sep != 0) %>%
    mutate(cumsum_sep = cumsum(replace_na(dec_sep, 0)))

  tolm_df <- tol_df %>%
    filter(dec_mar != 0) %>%
    mutate(cumsum_mar = cumsum(replace_na(dec_mar, 0)))

  tol_s <- tols_df$tol[which(tols_df$cumsum_sep >= max(tols_df$cumsum_sep) * k1)[1]]

  if (plot) {
    tols_plot <- ggplot(tols_df) +
      geom_point(aes(x = tol, y = cumsum_sep), color = 'grey25', alpha = 0.25) +
      geom_smooth(aes(x = tol, y = cumsum_sep),
                  method = "lm", formula = y ~ log(x), se = FALSE) +
      geom_vline(xintercept = tol_s, linetype = "dashed",
                 color = "red", linewidth = 1.2) +
      theme_bw() +
      theme(plot.title      = element_text(size = 20),
            axis.title.x    = element_text(size = 16),
            axis.text.x     = element_text(size = 14),
            axis.title.y    = element_text(size = 16),
            axis.text.y     = element_text(size = 14)) +
      labs(title = 'Latitude difference along fall equinox') +
      ylab('Latitude difference') + xlab('tol value')
  }

  if (length(tmar) >= 15) {
    tol_m <- tolm_df$tol[which(tolm_df$cumsum_mar >= max(tolm_df$cumsum_mar) * k2)[1]]
    if (plot) {
      tolm_plot <- ggplot(tolm_df) +
        geom_point(aes(x = tol, y = cumsum_mar), color = 'grey25', alpha = 0.25) +
        geom_smooth(aes(x = tol, y = cumsum_mar),
                    method = "lm", formula = y ~ log(x), se = FALSE) +
        geom_vline(xintercept = tol_m, linetype = "dashed",
                   color = "red", linewidth = 1.2) +
        theme_bw() +
        theme(plot.title   = element_text(size = 20),
              axis.title.x = element_text(size = 16),
              axis.text.x  = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y  = element_text(size = 14)) +
        labs(title = 'Latitude difference along spring equinox') +
        ylab('Latitude difference') + xlab('tol value')
    }
  } else {
    tol_m <- NA
    if (plot) {
      tolm_plot <- ggplot(tol_df) +
        labs(title = 'Latitude difference along spring equinox') +
        xlab('tol value') + ylab('Latitude difference') +
        xlim(c(0.00, 0.20)) + ylim(c(0, 100)) +
        theme_bw() +
        theme(plot.title   = element_text(size = 20),
              axis.title.x = element_text(size = 16),
              axis.text.x  = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y  = element_text(size = 14))
    }
  }

  cat('Tol values calculated\n')

  p <- if (plot) tols_plot / tolm_plot else NULL

  return(list(tol_sep = tol_s, tol_mar = tol_m, plot = p))
}


################################################################################

# Function to plot the SGAT track with equinox periods highlighted.
#
# @x0      : matrix of estimated positions (lon, lat columns)
# @time    : vector of twilight times corresponding to rows of x0
# @d       : number of days around equinox to highlight
# @Sep_eq  : Date of September equinox
# @Mar_eq  : Date of March equinox
# @world   : sf object for world map (from rnaturalearth)
# @col_lat : latitude of the colony
# @col_lon : longitude of the colony

plot_SGAT <- function(x0, time, d, Sep_eq, Mar_eq, col_lat, col_lon) {

  require(rnaturalearth)
  # Ensure x0 has named columns for ggplot
  colnames(x0) <- c("lon", "lat")

  require(ggplot2)

  xsep <- x0[which(as.Date(time) >= Sep_eq - d &
                     as.Date(time) <= Sep_eq + d), , drop = FALSE]
  xmar <- x0[which(as.Date(time) >= Mar_eq - d &
                     as.Date(time) <= Mar_eq + d), , drop = FALSE]

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  p <- ggplot() +
    geom_sf(data = world, fill = grey(0.9), color = grey(0.6), lwd = 0.2) +
    coord_sf(xlim = c(min(x0[, 1]), max(x0[, 1])),
             ylim = c(min(x0[, 2]), max(x0[, 2])), expand = TRUE) +
    geom_path(data  = data.frame(x0), aes(x = lon, y = lat), color = "firebrick") +
    geom_point(data = data.frame(x0),  aes(x = lon, y = lat),
               shape = 19, color = "firebrick") +
    geom_point(data = data.frame(xsep), aes(x = lon, y = lat),
               shape = 19, color = "cornflowerblue") +
    geom_point(data = data.frame(xmar), aes(x = lon, y = lat),
               shape = 19, color = "darkgreen") +
    geom_point(aes(x = col_lon, y = col_lat), shape = 16, size = 2.5) +
    xlab('Longitude') + ylab('Latitude') +
    theme_minimal() +
    theme(panel.grid   = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title   = element_text(size = 20),
          axis.title.x = element_text(size = 16),
          axis.text.x  = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y  = element_text(size = 14))

  return(p)
}
