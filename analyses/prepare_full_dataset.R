# =============================================================================
# PREPARE THE FULL SEABIRD DATASET (all species) FOR THE WORKFLOW
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Joins the GLS positions of all species used in the paper into
#              the two input files of the workflow, and runs basic checks.
#
#   Sources (all of them SGAT positions):
#     - Calonectris spp.: tracks listed in WAs_ids
#       (subset_CALtracks_WAclassified.RData), read from 'cal_dir'
#       (<geo_trip>_<year>_<Sp>_<Colony>_rawpoints[_cdays].csv)
#     - Hydrobates/*.csv         Hydrobates spp. (';' separator, ',' decimals)
#     - PUFPUF/*.csv             Puffinus puffinus
#     - PUFYEL/*.csv             Puffinus yelkouan (read only if the folder
#                                exists)
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

# Calonectris positions with colony anchoring ('_cdays' files: positions on
# colony days fixed at or near the colony) or without it ('_rawpoints' files).
# Main analysis: without anchoring (FALSE), as for the other species.
colony_anchoring <- FALSE

# Tracks excluded from the analysis (and reason)
exclude_tracks <- c('CJ636001',    # HYDCAS, Egg Island: not used in the paper
                    'BR928001001') # CALEDW, Montana Clara: not used in the paper

out_dir   <- 'data'
check_dir <- file.path(out_dir, 'checks')
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

# Species codes, latin and common names
species_names <- tribble(
  ~Sp,      ~Species,               ~Common_name,
  'CALBOR', 'Calonectris borealis', "Cory's shearwater",
  'CALDIO', 'Calonectris diomedea', "Scopoli's shearwater",
  'CALEDW', 'Calonectris edwardsii', 'Cape Verde shearwater',
  'PUFPUF', 'Puffinus puffinus',    'Manx shearwater',
  'PUFYEL', 'Puffinus yelkouan',    'Yelkouan shearwater',
  'HYDCAS', 'Hydrobates castro',    'Band-rumped storm petrel',
  'HYDJAB', 'Hydrobates jabejabe',  'Cape Verde storm petrel'
)


# =============================================================================
# 2. READ EACH SOURCE ----
# =============================================================================
# All sources are converted to the same columns:
#   geo_trip, Ring, Sp, Colony, Date_Time (UTC), lon, lat, source

# --- 2.1. Calonectris (SGAT csv) ---------------------------------------------
cal_env <- new.env()
load(file.path(raw_dir, 'subset_CALtracks_WAclassified.RData'), envir = cal_env)
cal_ids <- unique(as.character(cal_env$WAs_ids$geo_trip))

suffix <- ifelse(colony_anchoring, '_rawpoints_cdays.csv', '_rawpoints.csv')
cal_files <- list.files(cal_dir, pattern = paste0(suffix, '$'), full.names = TRUE)
cal_files <- cal_files[str_extract(basename(cal_files), '^[^_]+') %in% cal_ids]

missing_cal <- setdiff(cal_ids, str_extract(basename(cal_files), '^[^_]+'))
if (length(missing_cal) > 0) {
  warning('Calonectris tracks not found: ', paste(missing_cal, collapse = ', '))
}

read_sgat_csv <- function(f) {
  read_csv(f, show_col_types = FALSE,
           col_types = cols(Ring = 'c', geoID = 'c', geo_trip = 'c',
                            Date_Time = 'c')) %>%
    transmute(geo_trip = coalesce(geo_trip, geoID), Ring, Sp = Species,
              Colony, Date_Time = ymd_hms(Date_Time, tz = 'UTC'),
              lon, lat, source = basename(f))
}
cal <- map_dfr(cal_files, read_sgat_csv)

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

# --- 2.3. Puffinus puffinus and P. yelkouan (SGAT csv) -------------------
puf_files <- c(list.files(file.path(raw_dir, 'PUFPUF'), pattern = '\\.csv$',
                          full.names = TRUE),
               list.files(file.path(raw_dir, 'PUFYEL'), pattern = '\\.csv$',
                          full.names = TRUE))
puf <- map_dfr(puf_files, read_sgat_csv)

# --- 2.5. Join all -------------------------------------------------------------
all_pos <- bind_rows(cal, hyd, puf) %>%
  filter(!geo_trip %in% exclude_tracks) %>%
  mutate(Colony = recode(Colony, 'M.Clara' = 'MClara')) %>%
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
colonies_found <- all_pos %>% distinct(Sp, Colony) %>% arrange(Colony, Sp)

if (!file.exists(colonies_file)) {
  colonies_found %>%
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
# 4. METADATA (one row per track) ----
# =============================================================================

metadata <- all_pos %>%
  group_by(geo_trip) %>%
  summarise(Ring = first(Ring), Sp = first(Sp), Colony = first(Colony),
            n_positions = n(),
            start_date = min(Date_Time), end_date = max(Date_Time),
            duration_days = round(as.numeric(difftime(end_date, start_date,
                                                      units = 'days'))),
            source = first(source), .groups = 'drop') %>%
  left_join(species_names, by = 'Sp') %>%
  left_join(colonies %>% dplyr::select(Colony, Colony.lat, Colony.lon),
            by = 'Colony') %>%
  mutate(geo_alfanum = geo_trip) %>%
  relocate(geo_trip, Ring, Sp, Species, Common_name, Colony,
           Colony.lat, Colony.lon)


# =============================================================================
# 5. CHECKS ----
# =============================================================================

# --- 5.1. Distance to the colony at the start and end of each track ----------
dist_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  6371 * acos(pmin(1, sin(lat1 * rad) * sin(lat2 * rad) +
                     cos(lat1 * rad) * cos(lat2 * rad) * cos((lon2 - lon1) * rad)))
}

ends <- all_pos %>%
  group_by(geo_trip) %>%
  summarise(lon_first = median(head(lon, 10)), lat_first = median(head(lat, 10)),
            lon_last  = median(tail(lon, 10)), lat_last  = median(tail(lat, 10)),
            max_gap_days = round(max(as.numeric(diff(Date_Time), units = 'days')), 1),
            .groups = 'drop')

metadata <- metadata %>%
  left_join(ends, by = 'geo_trip') %>%
  mutate(start_dist_km = round(dist_km(lon_first, lat_first, Colony.lon, Colony.lat)),
         end_dist_km   = round(dist_km(lon_last,  lat_last,  Colony.lon, Colony.lat))) %>%
  dplyr::select(-lon_first, -lat_first, -lon_last, -lat_last)

# --- 5.2. Sample size by species and colony (compare with Table S1) ----------
sample_size <- metadata %>%
  count(Species, Common_name, Sp, Colony, name = 'n_tracks') %>%
  arrange(Species, Colony)
print(sample_size, n = Inf)
write_csv(sample_size, file.path(check_dir, 'sample_size.csv'))

# --- 5.3. Tracks to check -------------------------------------------------------
# Short tracks, long gaps, or tracks that do not start/end near the colony
# (median of the first/last 10 positions; GLS error is about 200 km).
to_check <- metadata %>%
  filter(duration_days < 240 | max_gap_days > 5 |
           start_dist_km > 1000 | end_dist_km > 1000 | is.na(Colony.lat)) %>%
  dplyr::select(geo_trip, Sp, Colony, n_positions, start_date, end_date,
                duration_days, max_gap_days, start_dist_km, end_dist_km)
print(to_check, n = Inf)
write_csv(to_check, file.path(check_dir, 'tracks_to_check.csv'))

# --- 5.4. Map of all tracks and colonies ---------------------------------------
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


# =============================================================================
# 6. SAVE ----
# =============================================================================

write_csv(all_pos %>% dplyr::select(geo_trip, Date_Time, lon, lat),
          file.path(out_dir, 'All_positions.csv'))
write_csv(metadata, file.path(out_dir, 'All_metadata.csv'))

cat('\nSaved', n_distinct(all_pos$geo_trip), 'tracks and', nrow(all_pos),
    'positions in', out_dir, '\n')
