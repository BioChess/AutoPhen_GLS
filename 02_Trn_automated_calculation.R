# =============================================================================
# TRN AUTOMATED CALCULATION
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Processes raw GLS light data into twilight events (.trn files).
#              Uses twilightCalc (GeoLight), twilight_cleanup and move_twilights
#              (SEATRACK) to calculate, clean and adjust twilight events.
#              Output .trn files are used as input for SGAT location modelling.
# =============================================================================


# =============================================================================
# 1. LOAD PACKAGES AND FUNCTIONS ----
# =============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  'readxl',                      # opening Excels and csv
  'tidyverse',                   # managing data
  'lubridate',                   # dealing with dates
  'evaluate',                    # check errors
  'GeoLight', 'BAStag', 'TwGeos', # geolocators
  'suncalc'                      # twilight times by location
)

# Detach purrr to avoid conflicts with base map()
if ("purrr" %in% .packages()) {
  detach("package:purrr", unload = TRUE)
}

# Load custom functions
source('functions/FUNCTION_twilight_cleanup.R')
source('functions/FUNCTION_move_twilights.R')
source('functions/01_FUNCTIONS_GLS_analysis.R')


# =============================================================================
# 2. READ METADATA AND CALIBRATION PARAMETERS ----
# =============================================================================

mdata     <- read_csv('input/GLS.files.csv',                    col_names = TRUE)
calibs    <- read_csv('output/Supervised_calibrations_SGAT.csv', col_names = TRUE)
calib.sum <- read_csv('output/Summary_SGAT_calibrations.csv',    col_names = TRUE)


# =============================================================================
# 3. MAIN LOOP: TWILIGHT CALCULATION PER GLS ----
# =============================================================================

gls.dir <- paste0(getwd(), '/input/GLS/')

for (i in seq_len(nrow(mdata))) {

  # --- 3.1. Read GLS information ---------------------------------------------

  info        <- mdata[i, ]
  lightfile   <- if_else(!is.na(info$File.adjlux), info$File.adjlux, info$File.lux)
  lightfile_c <- paste0(gls.dir, lightfile)

  if (!file.exists(lightfile_c) || is.na(lightfile_c)) {
    mdata[i, c(21:23)] <- list(NA, 'No data', 'No')
    message('Skipping geo (', i, '): light file not found')
    next
  }

  geoID    <- info$geo_alfanum
  geo_trip <- info$geo_trip
  Producer <- info$Producer
  col.name <- info$Colony
  col_lat  <- info$Colony.lat
  col_lon  <- info$Colony.lon

  start_date <- as.Date(info$start_date, tz = 'GMT')
  end_date   <- as.Date(info$end_date,   tz = 'GMT')
  year_rec   <- info$Year_rec

  cat('Processing geo (', i, ')', geo_trip, 'for year', year_rec, '\n')

  # Select zenith angle from calibration (or colony/year summary if uncalibrated)
  sun_angle <- select_sunAngle(
    producer  = Producer,
    geoID     = geoID,
    colony    = col.name,
    year      = year_rec,
    calib.sum = calib.sum,
    calibs    = calibs,
    param     = FALSE
  )

  # --- 3.2. Read light data --------------------------------------------------

  lu <- read_lightGLS(
    producer   = Producer,
    ID         = geoID,
    file       = lightfile_c,
    start = start_date,
    end   = end_date
  )

  if (is.null(lu) || nrow(lu) == 0) {
    mdata[i, c(21:23)] <- list(NA, 'No data', NA)
    message('Skipping geo (', i, '): no light data')
    next
  }

  # Recording interval in minutes
  li <- median(difftime(lu$dtime[2:nrow(lu)], lu$dtime[1:(nrow(lu) - 1)], units = 'mins'))

  # Start and end year of the recording
  start_year <- year(lu[1, 1])
  end_year   <- year(lu[nrow(lu), 1])

  # Skip GLS with less than one full year of data
  if (start_year == end_year) {
    mdata[i, c(21:23)] <- list(NA, 'Data', 'Remove')
    message('Skipping geo (', i, ') ', geo_trip, ': less than one year of data')
    next
  }

  # Equinox dates for this trip
  Sep_eq <- as.Date(paste0(start_year, '-09-21'))
  Mar_eq <- as.Date(paste0(end_year,   '-03-21'))

  # --- 3.3. Calculate twilight events ----------------------------------------

  twl <- try(
    twilightCalc(lu$dtime, lu$lux,
                 ask            = T,
                 preSelection   = TRUE,
                 LightThreshold = 2,
                 maxLight       = li),
    silent = TRUE
  )

  # For most recent R versions!!
  twl <- TwGeos::findTwilights(
    tagdata   = lu,
    threshold = 2,        
    include   = as.POSIXct(c(start_date, end_date), tz = "GMT"),
    extend    = 0,
    dark.min  = 15             # ignorar sombreados < 15 min
  )
  twl_gl <- export2GeoLight(twl)
  
  if (inherits(twl, 'try-error')) {
    mdata[i, c(21:23)] <- list(NA, 'Data. Problem with twilight calculation', NA)
    message('Skipping geo (', i, ') ', geo_trip, ': twilight calculation failed')
    next
  }

  # Fix type=1/2 misassignments (two consecutive same-type events)
  for (j in 1:(nrow(twl) - 2)) {
    if (sum(twl$type[j:(j + 2)]) == 3) twl$type[j + 1] <- 2
    if (sum(twl$type[j:(j + 2)]) == 6) twl$type[j + 1] <- 1
  }

  # GeoLight removes one extra minute vs BAStag/intiproc — compensate
  twl$tFirst[twl$type  == 2] <- twl$tFirst[twl$type  == 2] + 60
  twl$tSecond[twl$type == 1] <- twl$tSecond[twl$type == 1] + 60

  twl <- twl[, c(1:3)]

  # --- 3.4. Clean and adjust twilight events ---------------------------------

  par(mfrow = c(2, 1))

  twl2 <- twilight_cleanup(
    df                          = twl,
    breedingloc_lon             = col_lon,
    breedingloc_lat             = col_lat,
    months_breeding             = c(5, 6, 7, 8),
    months_extensive_legtucking = NA,
    show_plot                   = TRUE
  )
  twl2 <- distinct(twl2)

  twl3 <- move_twilights(
    df               = twl2,
    minutes_different = 15,
    sun              = sun_angle,
    show_plot        = TRUE
  )
  twl3 <- distinct(twl3)

  # Optional filters (uncomment to apply)
  twl3 <- twl3 %>%
    mutate(
      loess        = loessFilter(twl3, k = 3, plot = TRUE),
      speed_filter = distanceFilter(twl3, degElevation = sun_angle,
                                    distance = 500, units = "hour")
    )
  # %>% filter(loess == TRUE)
  # %>% filter(speed_filter == TRUE)

  # --- 3.5. Assign confidence values -----------------------------------------

  twl4 <- twl3 %>%
    mutate(
      Twilight   = tFirst,
      Rise       = ifelse(type == '1', TRUE, FALSE),
      confidence = 9,
      Julian     = yday(Twilight)
    )

  # Lower confidence around equinox periods
  twl4 <- twl4 %>%
    mutate(
      confidence = case_when(
        between(Twilight, Sep_eq - 20, Sep_eq + 20) ~ 2,
        between(Twilight, Sep_eq - 10, Sep_eq + 10) ~ 1,
        between(Twilight, Mar_eq - 20, Mar_eq + 20) ~ 2,
        between(Twilight, Mar_eq - 10, Mar_eq + 10) ~ 1,
        TRUE ~ confidence
      )
    )

  # --- 3.6. Save .trn files --------------------------------------------------

  trn <- twl4[, c(1:3, 6, 7, 9, 8, 4, 5)]

  file.trn <- paste0(geo_trip, '_', info$Sp, '_', year_rec, '_', info$Colony, '_automated.trn')

  if (!dir.exists('output/trn'))
    dir.create('output/trn', recursive = TRUE, showWarnings = FALSE)

  write.table(trn,
              file      = paste0('output/trn/', file.trn),
              sep       = ',',
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE)

  # Backup format compatible with Locator (BirdTracker) software
  trn.bak <- trn %>%
    mutate(
      date_time = format(Twilight, format = '%d/%m/%y %H:%M:%S'),
      type      = ifelse(Rise == TRUE, 'Sunrise', 'Sunset')
    ) %>%
    arrange(date_time) %>%
    dplyr::select(date_time, type, confidence)

  file.trn.bak <- paste0(geo_trip, '_', info$Sp, '_', year_rec, '_', info$Colony, '.trn')

  if (!dir.exists('output/trn_bak'))
    dir.create('output/trn_bak', recursive = TRUE, showWarnings = FALSE)

  write.table(trn.bak,
              file      = paste0('output/trn_bak/', file.trn.bak),
              sep       = ',',
              col.names = FALSE,
              row.names = FALSE,
              quote     = FALSE)

  cat('  Trn from geo (', i, ')', geo_trip, 'saved\n')

  # Update metadata with output file name
  mdata[i, c(23:25)] <- list(file.trn, 'Data', 'No')
}


# =============================================================================
# 4. SAVE UPDATED METADATA ----
# =============================================================================

write_csv(mdata, file = 'input/GLS.files2.csv')
cat('Updated metadata saved to input/GLS.files2.csv\n')
