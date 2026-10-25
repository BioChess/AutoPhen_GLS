# =============================================================================
# PUFFINUS YELKOUAN: SGAT POSITIONS FROM TWILIGHT FILES (.trn)
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: The Yelkouan shearwater tracks are only available as twilight
#              files (.trn, TransEdit format: date-time, Sunrise/Sunset,
#              confidence). This script estimates positions with the SGAT
#              threshold method, as in gls_preprocessing/03, but with a
#              standard zenith angle (no calibration data for these loggers).
#              Only the tracks from La Mola are used.
#
# Input:  <raw_dir>/PUFYEL/*_Mola.trn
#         data/colonies.csv (colony coordinates)
# Output: data/PUFYEL/<geo_trip>_PUFYEL_<Colony>_rawpoints.csv
#         (same columns as the SGAT '_rawpoints' files of Calonectris; read by
#         analyses/prepare_full_dataset.R)
# =============================================================================


# =============================================================================
# 1. PACKAGES, FUNCTIONS AND SETTINGS ----
# =============================================================================

# SGAT is only on GitHub. It is installed by 00_setup.R only when
# install_gls_preprocessing <- TRUE, so it is installed here if missing
# (same version as in 00_setup.R; on Windows it needs Rtools).
if (!requireNamespace('SGAT', quietly = TRUE)) {
  if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')
  remotes::install_github('SWotherspoon/SGAT@d0894847fb9b', upgrade = 'never')
}
library(SGAT)

pacman::p_load('tidyverse', 'lubridate', 'patchwork', 'ggplot2',
               'rnaturalearth', 'rnaturalearthdata')

# calc_tol() and plot_SGAT()
source('gls_preprocessing/functions/gls_preprocessing_functions.R')

raw_dir <- 'C:/Users/Diego Vicente/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/05_Algorithm_test/Seabird_tracks'
trn_pattern <- '_Mola\\.trn$'   # only the tracks from La Mola

# Colony name as written in data/colonies.csv
colony <- 'Mola'

# Standard zenith angle (no calibration available for these loggers).
# 97 is the mean zenith of the supervised calibrations of our Biotrack
# loggers. As a check, the script prints the zenith observed at the colony
# during the first 14 days after deployment.
zenith_standard <- 97

out_dir <- 'data/PUFYEL'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 2. COLONY COORDINATES ----
# =============================================================================

first_line <- readLines('data/colonies.csv', n = 1)
colonies <- if (grepl(';', first_line)) {
  read_csv2('data/colonies.csv', show_col_types = FALSE)
} else {
  read_csv('data/colonies.csv', show_col_types = FALSE)
}
col <- colonies %>% filter(Colony == colony)
if (nrow(col) != 1 || is.na(col$Colony.lat)) {
  stop('Colony "', colony, '" not found (or without coordinates) in data/colonies.csv')
}
col_lat <- col$Colony.lat
col_lon <- col$Colony.lon


# =============================================================================
# 3. POSITIONS FOR EACH TRN FILE ----
# =============================================================================

trn_files <- list.files(file.path(raw_dir, 'PUFYEL'), pattern = trn_pattern,
                        full.names = TRUE)
cat(length(trn_files), 'trn files\n')

for (f in trn_files) {

  # --- 3.1. Read twilights -----------------------------------------------------
  # File name: <geoID>_PUFYEL_<year>_<Colony>.trn; one trip per logger
  geoID    <- str_extract(basename(f), '^[^_]+')
  geo_trip <- paste0(geoID, '001')
  cat('\nTrack', geo_trip, '\n')

  twl <- read_csv(f, col_names = c('Twilight', 'Type', 'Confidence'),
                  col_types = cols(.default = 'c')) %>%
    mutate(Twilight = dmy_hms(Twilight, tz = 'UTC'),
           Rise     = Type == 'Sunrise') %>%
    filter(!is.na(Twilight)) %>%
    distinct(Twilight, .keep_all = TRUE) %>%
    arrange(Twilight)

  start_year <- year(min(twl$Twilight))
  end_year   <- year(max(twl$Twilight))

  # Check: zenith observed at the colony during the first 14 days
  first14 <- twl %>% filter(Twilight < min(Twilight) + days(14))
  z_col <- SGAT::zenith(SGAT::solar(first14$Twilight), col_lon, col_lat)
  cat('  Zenith at the colony (first 14 days): median', round(median(z_col), 1),
      '| standard zenith used:', zenith_standard, '\n')

  # --- 3.2. tol around the equinoxes (as in gls_preprocessing/03) --------------
  Sep_eq <- as.Date(paste0(start_year, '-09-21'))
  Mar_eq <- as.Date(paste0(end_year,   '-03-21'))
  tol_res <- calc_tol(twl, zenith = zenith_standard, n = 0.001, d = 20,
                      Sep_eq = Sep_eq, Mar_eq = Mar_eq,
                      k1 = 0.95, k2 = 0.5, plot = FALSE)

  # --- 3.3. SGAT threshold path ---------------------------------------------------
  if (!is.na(tol_res$tol_mar)) {
    # Different tol for each equinox: split the track on 1st January
    split_date <- as.POSIXct(paste0(end_year, '-01-01'), tz = 'UTC')
    twl_1  <- twl %>% filter(Twilight <= split_date)
    twl_2  <- twl %>% filter(Twilight >  split_date)
    path_1 <- SGAT::thresholdPath(twl_1$Twilight, twl_1$Rise,
                            zenith = zenith_standard, tol = tol_res$tol_sep)
    path_2 <- SGAT::thresholdPath(twl_2$Twilight, twl_2$Rise,
                            zenith = zenith_standard, tol = tol_res$tol_mar)
    x0   <- rbind(path_1$x, path_2$x)
    time <- c(path_1$time, path_2$time)
  } else {
    path <- SGAT::thresholdPath(twl$Twilight, twl$Rise,
                          zenith = zenith_standard, tol = tol_res$tol_sep)
    x0   <- path$x
    time <- path$time
  }

  # --- 3.4. Save (same columns as the Calonectris '_rawpoints' files) ------------
  # The first/last positions are set at the colony in prepare_full_dataset.R
  raw.points <- tibble(Ring = NA_character_, Species = 'PUFYEL',
                       Colony = colony, geoID = geoID, geo_trip = geo_trip,
                       Producer = NA_character_,
                       start_year = start_year, end_year = end_year,
                       zenith = zenith_standard,
                       tol_sep = tol_res$tol_sep, tol_mar = tol_res$tol_mar,
                       Date_Time = time,
                       lon = ((x0[, 1] + 180) %% 360) - 180,   # [-180, 180)
                       lat = x0[, 2])

  out_file <- file.path(out_dir, paste0(geo_trip, '_PUFYEL_', colony,
                                        '_rawpoints.csv'))
  write_csv(raw.points, out_file)
  cat('  ', nrow(raw.points), 'positions,', format(min(time), '%Y-%m-%d'), 'to',
      format(max(time), '%Y-%m-%d'), '| tol_sep =', tol_res$tol_sep,
      '| tol_mar =', tol_res$tol_mar, '\n')

  print(plot_SGAT(x0, time, d = 20, Sep_eq, Mar_eq, col_lat, col_lon) +
          ggtitle(geo_trip))
}
