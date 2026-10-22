# =============================================================================
# PREPARE GPS DATA FOR THE WORKFLOW
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Converts a GPS dataset in Seabird Tracking Database (STDB)
#              format into the two input files used by the workflow:
#                - positions: one row per position (geo_trip, Date_Time, lon, lat)
#                - metadata:  one row per track (ID, species, colony, colony
#                             coordinates)
#              GPS positions are resampled to a regular interval (by default,
#              two positions per day at 00:00 and 12:00 UTC, as in GLS data),
#              so that the same segmentation settings can be used for both
#              types of data and the computing time stays low.
#              If the data are already resampled, this step keeps them as
#              they are.
#
#              One track = one individual and one annual cycle. The track ID
#              (geo_trip) is built as <organismID>_<year of the first
#              position>, so a track joined from several deployments (or
#              datasets) is kept as a single track.
#
# Input:  input/ICHAUD_example.csv (STDB format)
# Output: input/GPS_positions_example.csv
#         input/GPS_metadata_example.csv
# =============================================================================


# =============================================================================
# 1. LOAD PACKAGES ----
# =============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load('tidyverse', 'lubridate')


# =============================================================================
# 2. SETTINGS ----
# =============================================================================

input_file     <- 'input/ICHAUD_example.csv'
out_positions  <- 'input/GPS_positions_example.csv'
out_metadata   <- 'input/GPS_metadata_example.csv'

# Resampling interval (hours). 12 h = two positions per day, as in GLS data.
resample_hours <- 12

# Maximum time difference (hours) between a target time (e.g. 12:00) and the
# closest GPS position. If there is no GPS position closer than this, the
# target time is left empty (no interpolation).
max_diff_hours <- 2


# =============================================================================
# 3. READ GPS DATA ----
# =============================================================================

gps <- read_csv(input_file, show_col_types = FALSE,
                col_types = cols(datetime = col_character())) %>%
  mutate(Date_Time = ymd_hms(datetime, tz = 'UTC')) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(Date_Time)) %>%
  arrange(organismID, Date_Time)

cat('GPS positions read:', nrow(gps), '\n',
    'Individuals:', n_distinct(gps$organismID), '\n',
    'From', format(min(gps$Date_Time)), 'to', format(max(gps$Date_Time)), '\n')


# =============================================================================
# 4. RESAMPLE TO A REGULAR INTERVAL ----
# =============================================================================
# For each individual we create a regular sequence of target times (00:00 and
# 12:00 UTC by default) and keep the GPS position closest to each target time.

resample_track <- function(df, every_h, max_h) {
  targets <- seq(floor_date(min(df$Date_Time), 'day'),
                 ceiling_date(max(df$Date_Time), 'day'),
                 by = paste(every_h, 'hours'))

  # Index of the closest GPS position to each target time
  idx <- findInterval(targets, df$Date_Time)
  idx_before <- pmax(idx, 1)
  idx_after  <- pmin(idx + 1, nrow(df))
  d_before <- abs(as.numeric(difftime(targets, df$Date_Time[idx_before], units = 'hours')))
  d_after  <- abs(as.numeric(difftime(df$Date_Time[idx_after], targets, units = 'hours')))
  closest  <- ifelse(d_before <= d_after, idx_before, idx_after)
  diff_h   <- pmin(d_before, d_after)

  df[closest, ] %>%
    mutate(target_time = targets, diff_h = diff_h) %>%
    filter(diff_h <= max_h)
}

gps_rs <- gps %>%
  group_split(organismID) %>%
  map(resample_track, every_h = resample_hours, max_h = max_diff_hours) %>%
  bind_rows()

cat('Positions after resampling:', nrow(gps_rs), '\n')


# =============================================================================
# 5. CREATE WORKFLOW INPUT FILES ----
# =============================================================================

# --- 5.1. Positions -----------------------------------------------------------
gps_rs <- gps_rs %>%
  group_by(organismID) %>%
  mutate(geo_trip = paste0(organismID, '_', year(min(Date_Time)))) %>%
  ungroup()

positions <- gps_rs %>%
  transmute(Ring      = organismID,
            Species   = speciesCode,
            Colony    = colonyName,
            geoID     = organismID,
            geo_trip  = geo_trip,
            Date_Time = Date_Time,
            lon       = longitude,
            lat       = latitude)

# --- 5.2. Metadata (one row per track) ----------------------------------------
metadata <- gps_rs %>%
  group_by(geo_trip) %>%
  summarise(Ring        = first(organismID),
            Colony      = first(colonyName),
            Colony.lat  = first(colonyLat),
            Colony.lon  = first(colonyLon),
            Sp          = first(speciesCode),
            geo_alfanum = first(organismID),
            Producer    = first(instrumentManufacturer),
            Device      = first(instrumentType),
            start_date  = min(Date_Time),
            end_date    = max(Date_Time),
            .groups = 'drop') %>%
  relocate(geo_trip, .after = end_date)

write_csv(positions, out_positions)
write_csv(metadata,  out_metadata)

cat('Saved:', out_positions, 'and', out_metadata, '\n')


# =============================================================================
# 6. QUICK CHECK: DISTANCE TO THE COLONY BY MONTH ----
# =============================================================================
# A migratory track should show months far from the colony (wintering).

dist_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  6371 * acos(pmin(1, sin(lat1 * rad) * sin(lat2 * rad) +
                     cos(lat1 * rad) * cos(lat2 * rad) * cos((lon2 - lon1) * rad)))
}

positions %>%
  left_join(metadata %>% dplyr::select(geo_trip, Colony.lon, Colony.lat),
            by = 'geo_trip') %>%
  mutate(month = floor_date(Date_Time, 'month'),
         dist  = dist_km(lon, lat, Colony.lon, Colony.lat)) %>%
  group_by(geo_trip, month) %>%
  summarise(n = n(),
            median_dist_km = round(median(dist), 1),
            max_dist_km    = round(max(dist), 1),
            .groups = 'drop') %>%
  print(n = Inf)
