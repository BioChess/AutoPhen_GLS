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
#     2. Calonectris: positions on colony days (identified from activity data
#        in the SGAT script) are removed. On these days the logger is in the
#        burrow and positions are wrong (often in Antarctica).
#     3. All species: positions outside a broad range of each species are
#        removed. This removes the same kind of errors in the species without
#        colony days information.
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

out_dir   <- 'data'
pufyel_dir <- file.path(out_dir, 'PUFYEL')
check_dir <- file.path(out_dir, 'checks')
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

# Colony days of Calonectris:
#   'remove' main analysis: positions on colony days are removed
#   'anchor' as in the submitted paper: positions on colony days are set at
#            the colony (confirmed days) or at 100 km (probable days)
#   'keep'   positions on colony days are kept (for sensitivity tests)
calonectris_colony_days <- 'remove'

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
  # Puffinus yelkouan (logger stopped on 16 March 2022)
  'BW101001001'
  # ADD incomplete tracks of the other species here
)

# Tracks excluded from the analysis (and reason)
exclude_tracks <- c('CJ636001',    # HYDCAS, Egg Island: not used in the paper
                    'BR928001001') # CALEDW, Montana Clara: not used in the paper

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
# '_rawpoints_cdays' file of the same track
read_cal <- function(f) {
  raw <- read_sgat_csv(f)
  f_cdays <- sub('_rawpoints\\.csv$', '_rawpoints_cdays.csv', f)
  if (!file.exists(f_cdays)) {
    warning('No _cdays file for ', basename(f), ': colony days not identified')
    return(mutate(raw, colony_day = FALSE))
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
         # 3. Outside the range of the species
         remove_range = !at_colony & !remove_cday &
           (lon < lon_min | lon > lon_max | lat < lat_min | lat > lat_max)) %>%
  ungroup()

cleaning <- all_pos %>%
  group_by(Sp) %>%
  summarise(n_tracks = n_distinct(geo_trip), n_positions = n(),
            removed_colony_days = sum(remove_cday),
            removed_out_of_range = sum(remove_range, na.rm = TRUE),
            .groups = 'drop') %>%
  mutate(percent_removed = round(100 * (removed_colony_days +
                                          removed_out_of_range) / n_positions, 1))
print(cleaning)
write_csv(cleaning, file.path(check_dir, 'cleaning_summary.csv'))

# Longest gap between positions before cleaning (gaps in the original data;
# removing colony days creates longer gaps during the breeding period)
gaps <- all_pos %>%
  group_by(geo_trip) %>%
  summarise(max_gap_days = round(max(as.numeric(diff(Date_Time), units = 'days')), 1),
            .groups = 'drop')

all_pos <- all_pos %>%
  filter(!remove_cday, !coalesce(remove_range, FALSE))


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
dist_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  6371 * acos(pmin(1, sin(lat1 * rad) * sin(lat2 * rad) +
                     cos(lat1 * rad) * cos(lat2 * rad) * cos((lon2 - lon1) * rad)))
}

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


# =============================================================================
# 7. SAVE ----
# =============================================================================

write_csv(all_pos %>% dplyr::select(geo_trip, Date_Time, lon, lat),
          file.path(out_dir, 'All_positions.csv'))
write_csv(metadata, file.path(out_dir, 'All_metadata.csv'))

cat('\nSaved', n_distinct(all_pos$geo_trip), 'tracks and', nrow(all_pos),
    'positions in', out_dir, '\n')
