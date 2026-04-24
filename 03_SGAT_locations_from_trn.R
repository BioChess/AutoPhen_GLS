# =============================================================================
# SGAT LOCATIONS FROM TRN FILES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Processes twilight events (.trn files) into geographic positions
#              using the SGAT threshold path method. Colony days information is
#              used to anchor known locations at the breeding site.
# =============================================================================


# =============================================================================
# 1. LOAD PACKAGES AND FUNCTIONS ----
# =============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  'readxl',                          # opening Excels and csv
  'tidyverse',                       # managing data
  'lubridate',                       # dealing with dates
  'SGAT',                            # geolocator threshold path model
  'evaluate',                        # check errors
  'geosphere',                       # destPoint for random colony offsets
  'rnaturalearth', 'rnaturalearthdata', # world map for plots
  'ggplot2', 'cowplot', 'patchwork'  # plots
)

# Detach purrr to avoid conflicts with base map()
if ("purrr" %in% .packages()) {
  detach("package:purrr", unload = TRUE)
}

source('functions/01_FUNCTIONS_GLS_analysis.R')


# =============================================================================
# 2. READ METADATA AND CALIBRATION PARAMETERS ----
# =============================================================================

mdata2    <- read_csv('input/GLS.files2.csv',                    col_names = TRUE)
calibs    <- read_csv('output/Supervised_calibrations_SGAT.csv', col_names = TRUE)
calib.sum <- read_csv('output/Summary_SGAT_calibrations.csv',    col_names = TRUE)
col.days  <- read_csv('input/colony_days.csv')


# =============================================================================
# 3. MAIN LOOP: SGAT POSITIONS PER GLS ----
# =============================================================================

raw_trips.list  <- vector('list', nrow(mdata2))
raw_trips2.list <- vector('list', nrow(mdata2))

for (i in seq_len(nrow(mdata2))) {

  # --- 3.1. Read GLS information ---------------------------------------------

  info     <- mdata2[i, ]
  trnfile  <- if_else(!is.na(info$File.trn), info$File.trn, NA_character_)

  if (is.na(trnfile)) {
    message('Skipping geo (', i, '): no trn file')
    next
  }

  trnfile_c <- paste0('output/', trnfile)

  if (!file.exists(trnfile_c)) {
    message('Skipping geo (', i, '): trn file not found at ', trnfile_c)
    next
  }

  ring      <- info$Ring
  sp        <- info$Sp
  geoID     <- info$geo_alfanum
  geo_trip  <- info$geo_trip
  Producer  <- info$Producer
  col.name  <- info$Colony
  col_lat   <- info$Colony.lat
  col_lon   <- info$Colony.lon
  year_rec  <- info$Year_rec

  start_date <- as.Date(info$start_date, tz = 'GMT')
  start_year <- year(start_date)
  end_date   <- as.Date(info$end_date, tz = 'GMT')
  end_year   <- year(end_date)

  cat('Processing geo (', i, ')', geo_trip, 'for year', year_rec,
      'from', col.name, '\n')

  # --- 3.2. Read trn file and calibration parameters ------------------------

  trn <- read_csv(trnfile_c, col_names = TRUE, show_col_types = FALSE)
  twl <- trn %>%
    filter(Twilight >= start_date & Twilight <= end_date) %>%
    distinct(Twilight, .keep_all = TRUE) %>%
    arrange(Twilight)

  if (nrow(twl) == 0) {
    message('Skipping geo (', i, ') ', geo_trip, ': no twilight events in date range')
    next
  }

  calib.param <- select_sunAngle(
    producer  = Producer,
    geoID     = geoID,
    colony    = col.name,
    year      = year_rec,
    calib.sum = calib.sum,
    calibs    = calibs,
    param     = TRUE
  )

  zenith  <- calib.param$zenith_mean
  zenith0 <- calib.param$zenith0_mean
  alpha   <- calib.param$alpha

  if (is.na(zenith)) {
    message('Skipping geo (', i, ') ', geo_trip, ': no zenith angle available')
    next
  }

  # --- 3.3. Determine tol parameter -----------------------------------------

  Sep_eq <- as.Date(paste0(start_year, '-09-21'))
  Mar_eq <- as.Date(paste0(end_year,   '-03-21'))

  tol_res <- tryCatch(
    calc_tol(twl, zenith = zenith, n = 0.001, d = 20,
             Sep_eq = Sep_eq, Mar_eq = Mar_eq,
             k1 = 0.95, k2 = 0.5, plot = FALSE),
    error = function(e) {
      message('  tol calculation failed: ', e$message)
      NULL
    }
  )

  if (is.null(tol_res)) next

  # --- 3.4. SGAT initial path -----------------------------------------------

  if (!is.null(tol_res$tol_mar) && !is.na(tol_res$tol_mar)) {
    # Split trip at Jan 1st to apply different tol for each equinox
    twl_1 <- twl %>% filter(Twilight >= start_date &
                               Twilight <= as.Date(paste0(end_year, '-01-01')))
    twl_2 <- twl %>% filter(Twilight >  as.Date(paste0(end_year, '-01-01')) &
                               Twilight <= end_date)

    path_1 <- thresholdPath(twl_1$Twilight, twl_1$Rise,
                             zenith = zenith, tol = tol_res$tol_sep)
    path_2 <- thresholdPath(twl_2$Twilight, twl_2$Rise,
                             zenith = zenith, tol = tol_res$tol_mar)

    x0   <- rbind(path_1$x, path_2$x)
    time <- c(path_1$time, path_2$time)

  } else {
    path <- thresholdPath(twl$Twilight, twl$Rise,
                          zenith = zenith, tol = tol_res$tol_sep)
    x0   <- path$x
    time <- path$time
  }

  colnames(x0) <- c('lon', 'lat')
  cat('  Initial path determined\n')

  # Plot raw track
  p1 <- plot_SGAT(x0, time, d = 20, Sep_eq, Mar_eq, col_lat, col_lon)
  print(p1)

  # --- 3.5. Save raw SGAT track ---------------------------------------------

  raw.points <- cbind(time, as.data.frame(x0))
  colnames(raw.points) <- c('Date_Time', 'lon', 'lat')

  raw.points <- bind_cols(
    Ring       = ring,       Species    = sp,
    Colony     = col.name,   geoID      = geoID,
    geo_trip   = geo_trip,   Producer   = Producer,
    start_year = start_year, end_year   = end_year,
    tol_sep    = tol_res$tol_sep,
    tol_mar    = tol_res$tol_mar,
    raw.points
  )

  mdata2[i, 'SGAT_raw'] <- 'Yes'

  # --- 3.6. Clean SGAT positions through colony days information ------------

  cdays <- col.days %>%
    dplyr::filter(geo_trip == .env$geo_trip) %>%
    filter(light == TRUE | wetdry == TRUE)

  auxi <- data.frame(date = as.Date(time)) %>%
    cbind(as.data.frame(x0), row = seq_len(nrow(x0))) %>%
    filter(date %in% cdays$date)

  auxi <- left_join(auxi, cdays[, c(3, 6, 9)], by = 'date')

  fixedx <- rep(0, nrow(x0))
  fixedx[auxi[is.na(auxi$wetdry),   ]$row] <- 3
  fixedx[auxi[auxi$wetdry == TRUE,  ]$row] <- 2
  fixedx[auxi[auxi$light  == TRUE,  ]$row] <- 1
  fixedx[1:2] <- 1  # first two positions always at colony

  # Last positions: colony if full deployment, free if resuscitated GLS
  if (end_date == as.Date(twl$Twilight[nrow(twl)])) {
    fixedx[(nrow(x0) - 1):nrow(x0)] <- 1
  } else {
    fixedx[(nrow(x0) - 1):nrow(x0)] <- 0
  }

  # Assign colony coordinates or nearby random point
  for (k in which(fixedx != 0)) {
    if (fixedx[k] == 1) {
      x0[k, 'lon'] <- col_lon
      x0[k, 'lat'] <- col_lat
    } else if (fixedx[k] == 2) {
      random_p     <- geosphere::destPoint(p = c(col_lon, col_lat),
                                           b = runif(1, 0, 360), d = 100 * 1000)
      x0[k, 'lon'] <- as.numeric(random_p[1])
      x0[k, 'lat'] <- as.numeric(random_p[2])
    }
  }

  # Plot colony-corrected track
  p2 <- plot_SGAT(x0, time, d = 20, Sep_eq, Mar_eq, col_lat, col_lon)
  print(p2)
  cat('  Colony days corrected\n')

  raw.points2 <- cbind(time, as.data.frame(x0))
  colnames(raw.points2) <- c('Date_Time', 'lon', 'lat')

  raw.points2 <- bind_cols(
    Ring       = ring,       Species    = sp,
    Colony     = col.name,   geoID      = geoID,
    geo_trip   = geo_trip,   Producer   = Producer,
    start_year = start_year, end_year   = end_year,
    tol_sep    = tol_res$tol_sep,
    tol_mar    = tol_res$tol_mar,
    raw.points2
  )

  raw_trips.list[[i]]  <- raw.points
  raw_trips2.list[[i]] <- raw.points2

  cat('Geo (', i, ')', geo_trip, 'done\n\n')
}


# =============================================================================
# 4. SAVE SGAT RAW TRACKS AND UPDATED METADATA ----
# =============================================================================

raw_trips  <- bind_rows(raw_trips.list)
raw_trips2 <- bind_rows(raw_trips2.list)

if (!dir.exists('output')) dir.create('output', recursive = TRUE)

write_csv(raw_trips,  file = 'output/Raw_SGAT_tracks_tol.csv')
write_csv(raw_trips2, file = 'output/Raw_SGAT_tracks_tol_cdays.csv')
write_csv(mdata2,     file = 'input/GLS.files3.csv')

cat('Raw SGAT tracks saved to output/\n')
cat('Updated metadata saved to input/GLS.files3.csv\n')
