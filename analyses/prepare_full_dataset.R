# =============================================================================
# PREPARE THE FULL SEABIRD DATASET (all species) FOR THE WORKFLOW
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Joins the GLS positions of all species used in the paper into
#              the two input files of the workflow, cleans them and runs basic
#              checks.
#
#   Sources (all of them SGAT positions without colony anchoring):
#     - Calonectris spp.: the 124 tracks used in the submitted paper (IDs
#       from CAL_subset_phenology_tracks_v9.5.csv), read from 'cal_dir':
#       <geo_trip>_<year>_<Sp>_<Colony>_rawpoints.csv        (positions)
#       <geo_trip>_<year>_<Sp>_<Colony>_rawpoints_cdays.csv  (to identify
#                                                             colony days)
#     - Hydrobates/*.csv         Hydrobates spp. (';' separator, ',' decimals)
#     - PUFPUF/*.csv             Puffinus puffinus
#     - data/PUFYEL/*.csv        Puffinus yelkouan (created from the .trn
#                                files with analyses/prepare_PUFYEL_positions.R)
#
#   Cleaning:
#     1. The first and last 'fix_ends_n' positions are set at the colony
#        (deployment and recovery; only the start in incomplete tracks).
#     2. Calonectris: positions on colony days (identified from light and
#        wet-dry data) are moved to a random point near the colony. On these
#        days the logger is in the burrow and positions are wrong (often in
#        Antarctica), but the bird is at the colony.
#     3. All species: positions outside a broad range of each species are
#        removed. This removes the same kind of errors in the species without
#        colony days information.
#     4. All species: spikes (short jumps far from the track that come back
#        after one to three positions) are removed.
#
#   Output (folder 'data/', NOT uploaded to GitHub):
#     - data/All_positions.csv   one row per position: geo_trip, Date_Time, lon, lat
#     - data/All_metadata.csv    one row per track: geo_trip, Ring, Sp,
#                                Species (latin), Common_name, Colony,
#                                Colony.lat, Colony.lon, start/end dates
#     - data/colonies.csv        colony coordinates (filled by the user)
#     - data/checks/             tables and map to check the data
# =============================================================================


# =============================================================================
# 1. PACKAGES AND SETTINGS ----
# =============================================================================

pacman::p_load('tidyverse', 'lubridate', 'sf', 'rnaturalearth',
               'rnaturalearthdata')

# Folders with the original files
raw_dir <- 'C:/Users/Diego Vicente/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/05_Algorithm_test/Seabird_tracks'
cal_dir <- 'C:/Users/Diego Vicente/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Results/SGAT/Raw_trips'
# File with the Calonectris tracks used in the submitted paper (only the IDs
# are used)
cal_ids_file <- 'C:/Users/Diego Vicente/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Results/Phenology/Algorithm/CAL_subset_phenology_tracks_v9.5.csv'
# Colony days of Calonectris (from light and wet-dry data). Only used for the
# tracks without a '_rawpoints_cdays.csv' file
colony_days_file <- 'C:/Users/Diego Vicente/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/Raw Data/colony_days_v3.csv'

out_dir   <- 'data'
pufyel_dir <- file.path(out_dir, 'PUFYEL')
check_dir <- file.path(out_dir, 'checks')
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

# Colony days of Calonectris (days at the colony from light and wet-dry data;
# GLS positions on these days are wrong because the logger is in the burrow):
#   'random' main analysis: positions on colony days are moved to a random
#            point within 'colony_radius_km' of the colony
#   'anchor' as in the submitted paper: positions on colony days are set at
#            the colony (confirmed days) or at exactly 100 km (probable days)
#   'remove' positions on colony days are removed
#   'keep'   positions on colony days are kept (for sensitivity tests)
calonectris_colony_days <- 'random'
colony_radius_km <- 100
set.seed(1)   # random points are the same every time the script is run

# Start and end of the tracks: the logger is deployed and recovered at the
# colony, so the first and last 'fix_ends_n' positions are set to the colony
# coordinates (the same rule used in the SGAT script). For incomplete tracks
# (e.g. logger failure before the return to the colony) only the start is
# fixed.
fix_ends_n <- 2
incomplete_tracks <- c(
  # Calonectris (from the '_cdays' files)
  '1605002002', '1644003001', '1644003002', '19084002001', 'BR50281001001',
  'BR50284001001', 'J459002001', 'X238001001', '1598002001', 'BR50345001001',
  '2280003001',
  # Puffinus puffinus (logger stopped a few days before the return)
  'BD732001001', '19205001001',
  # Puffinus yelkouan (logger stopped on 16 March 2022)
  'BW101001001'
  # ADD incomplete tracks of the other species here
)

# Tracks excluded from the analysis (and reason)
exclude_tracks <- c('CJ636001',    # HYDCAS, Egg Island: not used in the paper
                    'BR928001001') # CALEDW, Montana Clara: not used in the paper

# Spikes: positions more than 'spike_km' away from the running median of the
# 'spike_window' surrounding positions (odd number; 7 = 3.5 days with two
# positions per day) are removed (all species). NULL: no spike filter.
spike_km     <- 1000
spike_window <- 7

# Species codes, names and broad range (decimal degrees). Positions outside
# the range are removed. The range is broad on purpose (it includes the GLS
# error); it only removes positions that are clearly impossible.
species_info <- tribble(
  ~Sp,      ~Species,                ~Common_name,               ~lon_min, ~lon_max, ~lat_min, ~lat_max,
  'CALBOR', 'Calonectris borealis',  "Cory's shearwater",            -70,       60,      -50,       60,
  'CALDIO', 'Calonectris diomedea',  "Scopoli's shearwater",         -70,       60,      -50,       55,
  'CALEDW', 'Calonectris edwardsii', 'Cape Verde shearwater',        -70,        0,      -50,       40,
  'PUFPUF', 'Puffinus puffinus',     'Manx shearwater',              -75,       10,      -55,       70,
  'PUFYEL', 'Puffinus yelkouan',     'Yelkouan shearwater',          -10,       45,       28,       50,
  'HYDCAS', 'Hydrobates castro',     'Band-rumped storm petrel',    -100,        0,        0,       55,
  'HYDJAB', 'Hydrobates jabejabe',   'Cape Verde storm petrel',      -45,       10,      -40,       30
)


# =============================================================================
# 2. READ EACH SOURCE ----
# =============================================================================
# All sources are converted to the same columns:
#   geo_trip, Ring, Sp, Colony, Date_Time (UTC), lon, lat, colony_day, source

read_sgat_csv <- function(f) {
  read_csv(f, show_col_types = FALSE,
           col_types = cols(Ring = 'c', geoID = 'c', geo_trip = 'c',
                            Date_Time = 'c')) %>%
    transmute(geo_trip = coalesce(geo_trip, geoID), Ring, Sp = Species,
              Colony, Date_Time = ymd_hms(Date_Time, tz = 'UTC'),
              lon, lat, source = basename(f))
}

# --- 2.1. Calonectris (SGAT csv) ---------------------------------------------
cal_ids <- read_csv(cal_ids_file, col_select = 'ID', show_col_types = FALSE,
                    col_types = cols(ID = 'c')) %>%
  distinct(ID) %>% pull(ID)

cal_files <- list.files(cal_dir, pattern = '_rawpoints\\.csv$', full.names = TRUE)
cal_files <- cal_files[str_extract(basename(cal_files), '^[^_]+') %in% cal_ids]

missing_cal <- setdiff(cal_ids, str_extract(basename(cal_files), '^[^_]+'))
if (length(missing_cal) > 0) {
  warning('Calonectris tracks not found: ', paste(missing_cal, collapse = ', '))
}

# Colony days: positions that differ between the '_rawpoints' and the
# '_rawpoints_cdays' file of the same track. If there is no '_cdays' file,
# the colony days are taken from 'colony_days_file' (days with light or
# wet-dry evidence of being at the colony, the same rule as the SGAT script).
colony_days <- read_csv(colony_days_file, show_col_types = FALSE,
                        col_types = cols(geoID = 'c', geo_trip = 'c',
                                         date = 'D')) %>%
  filter(light %in% TRUE | wetdry %in% TRUE) %>%
  distinct(geo_trip, date)

read_cal <- function(f) {
  raw <- read_sgat_csv(f)
  f_cdays <- sub('_rawpoints\\.csv$', '_rawpoints_cdays.csv', f)
  if (!file.exists(f_cdays)) {
    id <- raw$geo_trip[1]
    cday_dates <- colony_days$date[colony_days$geo_trip == id]
    if (length(cday_dates) == 0) {
      warning('No colony days found for ', id)
    } else {
      message('No _cdays file for ', id, ': colony days from ',
              basename(colony_days_file))
    }
    if (calonectris_colony_days == 'anchor') {
      warning('Positions of ', id, ' not anchored (no _cdays file)')
    }
    return(mutate(raw, colony_day = as.Date(Date_Time) %in% cday_dates))
  }
  cdays <- read_sgat_csv(f_cdays)
  if (nrow(cdays) != nrow(raw) || any(cdays$Date_Time != raw$Date_Time)) {
    stop('Different positions in ', basename(f), ' and its _cdays file')
  }
  raw <- raw %>%
    mutate(colony_day = abs(lon - cdays$lon) > 1e-6 | abs(lat - cdays$lat) > 1e-6)
  if (calonectris_colony_days == 'anchor') {
    raw <- raw %>% mutate(lon = cdays$lon, lat = cdays$lat)
  }
  raw
}
cal <- map_dfr(cal_files, read_cal)

# --- 2.2. Hydrobates (SGAT csv; ';' separator and ',' decimals) ---------------
hyd_files <- list.files(file.path(raw_dir, 'Hydrobates'), pattern = '\\.csv$',
                        full.names = TRUE)
hyd <- map_dfr(hyd_files, function(f) {
  read_delim(f, delim = ';', locale = locale(decimal_mark = ','),
             show_col_types = FALSE,
             col_types = cols(Ring = 'c', geoID = 'c', geo_trip = 'c',
                              Date_Time = 'c')) %>%
    transmute(geo_trip = geoID,              # geo_trip is empty in these files
              Ring, Sp = Species, Colony,
              Date_Time = ymd_hms(Date_Time, tz = 'UTC'),
              lon, lat, source = basename(f))
})

# --- 2.3. Puffinus puffinus and P. yelkouan (SGAT csv) -------------------------
puf_files <- c(list.files(file.path(raw_dir, 'PUFPUF'), pattern = '\\.csv$',
                          full.names = TRUE),
               list.files(pufyel_dir, pattern = '_rawpoints\\.csv$',
                          full.names = TRUE))
if (length(list.files(pufyel_dir)) == 0) {
  warning('No P. yelkouan files: run analyses/prepare_PUFYEL_positions.R first')
}
puf <- map_dfr(puf_files, read_sgat_csv)

# --- 2.4. Join all -------------------------------------------------------------
all_pos <- bind_rows(cal, hyd, puf) %>%
  mutate(colony_day = coalesce(colony_day, FALSE)) %>%
  filter(!geo_trip %in% exclude_tracks) %>%
  mutate(Colony = recode(Colony, 'M.Clara' = 'MClara'),
         lon = ((lon + 180) %% 360) - 180) %>%   # longitudes in [-180, 180)
  filter(!is.na(lon), !is.na(lat), !is.na(Date_Time)) %>%
  distinct(geo_trip, Date_Time, .keep_all = TRUE) %>%
  arrange(Sp, geo_trip, Date_Time)


# =============================================================================
# 3. COLONY COORDINATES ----
# =============================================================================
# data/colonies.csv must contain one row per colony with its coordinates.
# The first time, a template is created: please fill Colony.lat and Colony.lon
# (decimal degrees, WGS84) and run this script again.

colonies_file <- file.path(out_dir, 'colonies.csv')

if (!file.exists(colonies_file)) {
  all_pos %>%
    distinct(Colony) %>%
    mutate(Colony.lat = NA_real_, Colony.lon = NA_real_, Notes = '') %>%
    write_csv(colonies_file)
  message('Template created: ', colonies_file,
          '. Fill the coordinates and run the script again.')
}
# The file can be saved with ',' or ';' as separator (e.g. Excel in Spanish)
first_line <- readLines(colonies_file, n = 1)
colonies <- if (grepl(';', first_line)) {
  read_csv2(colonies_file, show_col_types = FALSE)
} else {
  read_csv(colonies_file, show_col_types = FALSE)
}

missing_col <- setdiff(unique(all_pos$Colony),
                       colonies$Colony[!is.na(colonies$Colony.lat)])
if (length(missing_col) > 0) {
  warning('Colonies without coordinates: ', paste(missing_col, collapse = ', '))
}


# =============================================================================
# 4. CLEANING ----
# =============================================================================

# Random point within 'radius_km' of a point (uniform over the area of the
# circle). Used for the positions on colony days.
random_point <- function(lon, lat, radius_km) {
  n     <- length(lon)
  d     <- radius_km * sqrt(runif(n)) / 6371          # angular distance
  b     <- runif(n, 0, 2 * pi)                       # bearing
  rad   <- pi / 180
  lat1  <- lat * rad; lon1 <- lon * rad
  lat2  <- asin(sin(lat1) * cos(d) + cos(lat1) * sin(d) * cos(b))
  lon2  <- lon1 + atan2(sin(b) * sin(d) * cos(lat1),
                        cos(d) - sin(lat1) * sin(lat2))
  list(lon = lon2 / rad, lat = lat2 / rad)
}

# Great-circle distance (km)
dist_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  6371 * acos(pmin(1, sin(lat1 * rad) * sin(lat2 * rad) +
                     cos(lat1 * rad) * cos(lat2 * rad) * cos((lon2 - lon1) * rad)))
}

# Spikes: positions that are more than 'spike_km' away from the running
# median of the surrounding positions ('spike_window' positions centred on
# each one). The running median follows the real movement of the bird (also
# during fast migration), but not short jumps of one to three positions that
# come back. Positions fixed at the colony are never removed.
find_spikes <- function(lon, lat, fixed, spike_km, spike_window) {
  if (length(lon) < spike_window) return(rep(FALSE, length(lon)))
  lon_med <- runmed(lon, k = spike_window, endrule = 'median')
  lat_med <- runmed(lat, k = spike_window, endrule = 'median')
  dist_km(lon, lat, lon_med, lat_med) > spike_km & !fixed
}

all_pos <- all_pos %>%
  left_join(colonies %>% dplyr::select(Colony, Colony.lat, Colony.lon),
            by = 'Colony') %>%
  left_join(species_info %>% dplyr::select(Sp, lon_min, lon_max, lat_min, lat_max),
            by = 'Sp') %>%
  group_by(geo_trip) %>%
  mutate(k = row_number(), n = n(),
         # 1. First and last positions at the colony
         at_colony = k <= fix_ends_n |
           (k > n - fix_ends_n & !geo_trip %in% incomplete_tracks),
         lon = ifelse(at_colony & !is.na(Colony.lon), Colony.lon, lon),
         lat = ifelse(at_colony & !is.na(Colony.lat), Colony.lat, lat),
         # 2. Colony days (Calonectris)
         remove_cday = colony_day & !at_colony &
           calonectris_colony_days == 'remove',
         move_cday = colony_day & !at_colony &
           calonectris_colony_days == 'random') %>%
  ungroup()

# 2b. Colony days moved to a random point near the colony (option 'random')
if (any(all_pos$move_cday)) {
  rp <- with(all_pos[all_pos$move_cday, ],
             random_point(Colony.lon, Colony.lat, colony_radius_km))
  all_pos$lon[all_pos$move_cday] <- rp$lon
  all_pos$lat[all_pos$move_cday] <- rp$lat
}

# 3. Outside the range of the species (positions at or near the colony are
#    never removed)
all_pos <- all_pos %>%
  mutate(remove_range = !at_colony & !move_cday & !remove_cday &
           (lon < lon_min | lon > lon_max | lat < lat_min | lat > lat_max))

# 4. Spikes (checked on the positions that remain after steps 2 and 3)
all_pos <- all_pos %>%
  mutate(remove_spike = FALSE)
if (!is.null(spike_km)) {
  keep <- !all_pos$remove_cday & !coalesce(all_pos$remove_range, FALSE)
  all_pos$remove_spike[keep] <- all_pos[keep, ] %>%
    group_by(geo_trip) %>%
    mutate(spike = find_spikes(lon, lat, at_colony | move_cday, spike_km,
                               spike_window)) %>%
    pull(spike)
}

cleaning <- all_pos %>%
  group_by(Sp) %>%
  summarise(n_tracks = n_distinct(geo_trip), n_positions = n(),
            moved_colony_days = sum(move_cday),
            removed_colony_days = sum(remove_cday),
            removed_out_of_range = sum(remove_range, na.rm = TRUE),
            removed_spikes = sum(remove_spike),
            .groups = 'drop') %>%
  mutate(percent_removed = round(100 * (removed_colony_days +
                                          removed_out_of_range +
                                          removed_spikes) / n_positions, 1))
print(cleaning)
write_csv(cleaning, file.path(check_dir, 'cleaning_summary.csv'))

# Longest gap between positions before cleaning (gaps in the original data;
# removing colony days creates longer gaps during the breeding period)
gaps <- all_pos %>%
  group_by(geo_trip) %>%
  summarise(max_gap_days = round(max(as.numeric(diff(Date_Time), units = 'days')), 1),
            .groups = 'drop')

pos_before_cleaning <- all_pos   # kept for the plots of each track (6.5)

all_pos <- all_pos %>%
  filter(!remove_cday, !coalesce(remove_range, FALSE), !remove_spike)


# =============================================================================
# 5. METADATA (one row per track) ----
# =============================================================================

metadata <- all_pos %>%
  group_by(geo_trip) %>%
  summarise(Ring = first(Ring), Sp = first(Sp), Colony = first(Colony),
            n_positions = n(),
            start_date = min(Date_Time), end_date = max(Date_Time),
            duration_days = round(as.numeric(difftime(end_date, start_date,
                                                      units = 'days'))),
            source = first(source), .groups = 'drop') %>%
  left_join(species_info %>% dplyr::select(Sp, Species, Common_name), by = 'Sp') %>%
  left_join(colonies %>% dplyr::select(Colony, Colony.lat, Colony.lon),
            by = 'Colony') %>%
  mutate(geo_alfanum = geo_trip,
         incomplete = geo_trip %in% incomplete_tracks) %>%
  relocate(geo_trip, Ring, Sp, Species, Common_name, Colony,
           Colony.lat, Colony.lon)


# =============================================================================
# 6. CHECKS ----
# =============================================================================

# --- 6.1. Distance to the colony at the start and end of each track ----------

# Median of the first/last 10 positions that were not fixed at the colony:
# large distances suggest incomplete tracks
ends <- all_pos %>%
  filter(!at_colony) %>%
  group_by(geo_trip) %>%
  summarise(lon_first = median(head(lon, 10)), lat_first = median(head(lat, 10)),
            lon_last  = median(tail(lon, 10)), lat_last  = median(tail(lat, 10)),
            .groups = 'drop')
metadata <- metadata %>%
  left_join(ends, by = 'geo_trip') %>%
  left_join(gaps, by = 'geo_trip') %>%
  mutate(start_dist_km = round(dist_km(lon_first, lat_first, Colony.lon, Colony.lat)),
         end_dist_km   = round(dist_km(lon_last,  lat_last,  Colony.lon, Colony.lat))) %>%
  dplyr::select(-lon_first, -lat_first, -lon_last, -lat_last)

# --- 6.2. Sample size by species and colony (compare with Table S1) ----------
sample_size <- metadata %>%
  count(Species, Common_name, Sp, Colony, name = 'n_tracks') %>%
  arrange(Species, Colony)
print(sample_size, n = Inf)
write_csv(sample_size, file.path(check_dir, 'sample_size.csv'))

# --- 6.3. Tracks to check -------------------------------------------------------
# Short tracks, long gaps in the original data, or tracks that do not start/end near the colony
# (GLS error is about 200 km, but birds can forage far from the colony).
to_check <- metadata %>%
  filter(duration_days < 240 | max_gap_days > 5 |
           start_dist_km > 1000 | end_dist_km > 1000 | is.na(Colony.lat)) %>%
  dplyr::select(geo_trip, Sp, Colony, incomplete, n_positions, start_date,
                end_date, duration_days, max_gap_days, start_dist_km,
                end_dist_km)
print(to_check, n = Inf)
write_csv(to_check, file.path(check_dir, 'tracks_to_check.csv'))

# --- 6.4. Map of all tracks and colonies ---------------------------------------
world <- ne_countries(scale = 'medium', returnclass = 'sf')
map_all <- ggplot() +
  geom_sf(data = world, fill = 'grey85', colour = 'grey60', linewidth = 0.1) +
  geom_path(data = all_pos, aes(lon, lat, group = geo_trip, colour = Sp),
            linewidth = 0.2, alpha = 0.5) +
  geom_point(data = filter(colonies, !is.na(Colony.lat)),
             aes(Colony.lon, Colony.lat), shape = 23, size = 3,
             fill = 'orchid3') +
  coord_sf(xlim = range(all_pos$lon), ylim = range(all_pos$lat)) +
  facet_wrap(~ Sp) +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_bw() + theme(legend.position = 'none')
ggsave(file.path(check_dir, 'map_all_tracks.png'), map_all,
       width = 30, height = 20, units = 'cm', dpi = 150)


# --- 6.5. Plot of each track ------------------------------------------------------
# Map, longitude and latitude of each track, with the positions removed in the
# cleaning. Only for the species in 'plot_species' (NULL: no plots).
plot_species <- c('PUFPUF', 'PUFYEL', 'HYDCAS', 'HYDJAB')

plot_track <- function(id) {
  d <- pos_before_cleaning %>%
    filter(geo_trip == id) %>%
    mutate(status = case_when(at_colony    ~ 'Fixed at colony',
                              move_cday    ~ 'Colony day (moved)',
                              remove_cday  ~ 'Removed: colony day',
                              coalesce(remove_range, FALSE) ~ 'Removed: out of range',
                              remove_spike ~ 'Removed: spike',
                              TRUE         ~ 'Kept'),
           status = factor(status, levels = c('Kept', 'Fixed at colony',
                                       'Colony day (moved)',
                                              'Removed: colony day',
                                              'Removed: out of range',
                                              'Removed: spike')))
  kept <- filter(d, status %in% c('Kept', 'Fixed at colony', 'Colony day (moved)'))
  info <- filter(metadata, geo_trip == id)
  cols <- c('Kept' = 'grey20', 'Fixed at colony' = 'orchid3',
            'Colony day (moved)' = 'darkgreen',
            'Removed: colony day' = 'dodgerblue', 'Removed: out of range' = 'red',
            'Removed: spike' = 'orange')
  lims_x <- range(kept$lon) + c(-5, 5); lims_y <- range(kept$lat) + c(-5, 5)

  p_map <- ggplot() +
    geom_sf(data = world, fill = 'grey85', colour = 'grey60', linewidth = 0.1) +
    geom_path(data = kept, aes(lon, lat), colour = 'grey50', linewidth = 0.3) +
    geom_point(data = d, aes(lon, lat, colour = status), size = 0.8) +
    geom_point(data = info, aes(Colony.lon, Colony.lat), shape = 23, size = 3,
               fill = 'orchid3') +
    coord_sf(xlim = lims_x, ylim = lims_y) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = NULL, colour = NULL) +
    theme_bw() + theme(legend.position = 'bottom')

  d_long <- d %>%
    dplyr::select(Date_Time, lon, lat, status) %>%
    pivot_longer(c(lon, lat), names_to = 'coord', values_to = 'value')
  col_ref <- tibble(coord = c('lon', 'lat'),
                    value = c(info$Colony.lon, info$Colony.lat))
  p_time <- ggplot(d_long, aes(Date_Time, value)) +
    geom_hline(data = col_ref, aes(yintercept = value), colour = 'orchid3',
               linetype = 'dashed') +
    geom_line(data = filter(d_long, status %in% c('Kept', 'Fixed at colony',
                                                  'Colony day (moved)')),
              colour = 'grey50', linewidth = 0.3) +
    geom_point(aes(colour = status), size = 0.6) +
    facet_wrap(~ coord, ncol = 1, scales = 'free_y') +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw() + theme(legend.position = 'none')

  (p_map | p_time) +
    plot_layout(widths = c(1, 1.3)) +
    plot_annotation(title = paste0(id, ' - ', info$Sp, ' - ', info$Colony,
                                   ifelse(info$incomplete, ' (incomplete)', '')),
                    subtitle = paste0(format(info$start_date, '%Y-%m-%d'), ' to ',
                                      format(info$end_date, '%Y-%m-%d'),
                                      ' | start/end distance to colony: ',
                                      info$start_dist_km, ' / ', info$end_dist_km,
                                      ' km'))
}

if (length(plot_species) > 0) {
  pacman::p_load('patchwork')
  for (sp in plot_species) {
    sp_dir <- file.path(check_dir, 'tracks', sp)
    dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)
    for (id in metadata$geo_trip[metadata$Sp == sp]) {
      ggsave(file.path(sp_dir, paste0(id, '.png')), plot_track(id),
             width = 30, height = 15, units = 'cm', dpi = 120)
    }
  }
}


# =============================================================================
# 7. SAVE ----
# =============================================================================

write_csv(all_pos %>% dplyr::select(geo_trip, Date_Time, lon, lat),
          file.path(out_dir, 'All_positions.csv'))
write_csv(metadata, file.path(out_dir, 'All_metadata.csv'))

cat('\nSaved', n_distinct(all_pos$geo_trip), 'tracks and', nrow(all_pos),
    'positions in', out_dir, '\n')
