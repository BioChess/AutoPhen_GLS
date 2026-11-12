# =============================================================================
# SIMULATIONS 01: SIMULATED MIGRATORY TRACKS AT REAL SCALES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Simulates tracks of a seabird that breeds in a colony
#              (Veneguera, Gran Canaria), migrates to one or two wintering
#              areas and comes back. Movement is simulated with marcher
#              (Gurarie et al. 2017) around a moving centre of attraction, in
#              km and days. Observation error of different tracking methods is
#              added to the true positions.
#
#   Scenarios (full factorial, n_rep tracks each):
#     distance  short (800-2000 km), medium (3000-5000 km), long (6000-9000 km)
#     movement  base      one breeding area and one wintering area
#               foraging  central-place foraging trips (4-20 days) from the colony during
#                         the breeding periods
#               stopover  a stopover during the outward migration
#               two_W     two wintering areas
#     error     none      no observation error (e.g. GPS)
#               gls       GLS-like error (larger in latitude, autocorrelated,
#                         latitude error inflated around the equinoxes,
#                         missing positions); migrations away from the
#                         equinoxes
#               gls_eq    the same error, with migrations overlapping the
#                         equinoxes
#
# Output (simulations/data/):
#   sim_positions.csv        observed positions (geo_trip, Date_Time, lon, lat):
#                            the same format as data/All_positions.csv
#   sim_metadata.csv         one row per track (format of data/All_metadata.csv
#                            + scenario)
#   sim_truth_dates.csv      true dates and parameters of each track
#   sim_truth_positions.csv  true position and true class (B, mig, W) at each
#                            observed time
#   checks/                  maps of example tracks and summary of the error
#
# Run from the project root. Needs the package marcher (CRAN).
# =============================================================================


# =============================================================================
# 1. PACKAGES AND SETTINGS ----
# =============================================================================

pacman::p_load('tidyverse', 'sf', 'marcher', 'patchwork',
               'rnaturalearth', 'rnaturalearthdata')
source('simulations/simulation_functions.R')

out_dir <- file.path('simulations', 'data')
dir.create(file.path(out_dir, 'checks'), recursive = TRUE, showWarnings = FALSE)

p <- list(
  # Colony (Veneguera, Gran Canaria)
  colony = 'Veneguera', col_lon = -15.789, col_lat = 27.844,
  positions_per_day = 2,
  start_date = '2021-07-01',          # first position (breeding period)

  # Number of tracks per scenario
  n_rep = 25,

  # --- Movement --------------------------------------------------------------------
  # Distance from the colony to the (first) wintering area, by class (km)
  distance_km = list(short = c(800, 2000), medium = c(3000, 5000),
                     long = c(6000, 9000)),
  bearing_deg = c(190, 230),          # direction of the wintering area (from north;
                                      # south-west, along the Atlantic)
  speed_km_day = c(300, 600),         # migration speed
  range_radius_km = c(100, 250),      # radius with 95% of the positions around
                                      # the centre (marcher A = pi * r^2)
  tau_z_days = c(2.5, 3.5),           # position autocorrelation (marcher tau.z)
  tau_v_days = c(0.02, 0.03),         # velocity autocorrelation (marcher tau.v)
  min_wintering_days = 30,            # minimum stay in each wintering area
  after_arrival_days = c(15, 45),     # breeding period after the return
  # foraging trips (movement 'foraging')
  trip_gap_days = c(1, 15),           # time at the colony between trips
  trip_days = c(4, 20),               # duration of a trip
  trip_km = c(150, 1000),             # maximum distance from the colony
  # stopover (movement 'stopover')
  stopover_days = c(5, 15),
  stopover_position = c(0.3, 0.7),    # proportion of the route
  # second wintering area (movement 'two_W')
  w2_km = c(500, 1500),               # distance between the two wintering areas

  # --- Timing (month-day ranges) ----------------------------------------------------
  # normal:  migrations away from the equinoxes (+-20 days)
  # equinox: outward and return migrations overlap the equinoxes
  timing = list(normal  = list(B1_dep = c('10-20', '11-15'), B2_arr = c('02-05', '02-25')),
                equinox = list(B1_dep = c('09-10', '10-02'), B2_arr = c('03-08', '03-28'))),

  # --- Observation error ----------------------------------------------------------
  error = list(
    # GLS: mean error of about 150-200 km, larger in latitude (Phillips et al.
    # 2004; Halpin et al. 2021); errors correlated between consecutive positions;
    # latitude almost uninformative within +-eq_days of each equinox
    gls = list(sd_lon_km = 90, sd_lat_km = 180, phi = 0.5,
               eq_days = 15, eq_lat_factor = 4, eq_phi = 0.9,
               p_missing = 0.03, gaps_per_track = 1, gap_days = c(1, 3)),
    # Argos PTT (optional block): error of a few km
    ptt = list(sd_km = 5)),
  include_ptt = FALSE
)


# =============================================================================
# 2. SCENARIOS ----
# =============================================================================

scenarios <- expand_grid(distance = c('short', 'medium', 'long'),
                         movement = c('base', 'foraging', 'stopover', 'two_W'),
                         error = c('none', 'gls', 'gls_eq'),
                         rep = seq_len(p$n_rep))
if (p$include_ptt) {
  scenarios <- bind_rows(scenarios,
                         expand_grid(distance = c('short', 'medium', 'long'),
                                     movement = c('base', 'foraging', 'stopover', 'two_W'),
                                     error = 'ptt', rep = 1:5))
}
scenarios <- scenarios %>%
  mutate(n = row_number(),
         sim_id = sprintf('SIM%04d', n),
         seed = n)
cat(nrow(scenarios), 'tracks to simulate\n')


# =============================================================================
# 3. SIMULATION ----
# =============================================================================

positions <- list(); truth_pos <- list(); truth_dates <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  s  <- simulate_track(sc, p)
  obs <- add_error(s$true, sc$error, p, sc$seed)

  # The first and last two positions are set at the colony (known deployment
  # and recovery), as in the real data (analyses/prepare_full_dataset.R)
  n <- nrow(obs)
  obs$x[c(1:2, (n - 1):n)] <- 0; obs$y[c(1:2, (n - 1):n)] <- 0

  ll  <- km_to_lonlat(obs$x, obs$y, p$col_lon, p$col_lat)
  llt <- km_to_lonlat(obs$x_true, obs$y_true, p$col_lon, p$col_lat)
  positions[[i]] <- tibble(geo_trip = sc$sim_id, Date_Time = obs$Date_Time,
                           lon = ll$lon, lat = ll$lat)
  truth_pos[[i]] <- tibble(geo_trip = sc$sim_id, Date_Time = obs$Date_Time,
                           lon_true = llt$lon, lat_true = llt$lat,
                           error_km = sqrt((obs$x - obs$x_true)^2 + (obs$y - obs$y_true)^2),
                           stage = obs$stage, leg = obs$leg)
  truth_dates[[i]] <- s$dates
  if (i %% 50 == 0) cat('  ', i, 'tracks\n')
}
positions   <- bind_rows(positions)
truth_pos   <- bind_rows(truth_pos)
truth_dates <- bind_rows(truth_dates)

metadata <- positions %>%
  group_by(geo_trip) %>%
  summarise(n_positions = n(), start_date = min(Date_Time), end_date = max(Date_Time),
            .groups = 'drop') %>%
  mutate(Ring = geo_trip, Sp = 'SIM', Species = 'Simulated',
         Common_name = 'Simulated seabird', Colony = p$colony,
         Colony.lat = p$col_lat, Colony.lon = p$col_lon,
         duration_days = round(as.numeric(difftime(end_date, start_date, units = 'days')))) %>%
  left_join(truth_dates %>% select(geo_trip = sim_id, distance, movement, error),
            by = 'geo_trip') %>%
  select(geo_trip, Ring, Sp, Species, Common_name, Colony, Colony.lat, Colony.lon,
         n_positions, start_date, end_date, duration_days, distance, movement, error)

write_csv(positions,   file.path(out_dir, 'sim_positions.csv'))
write_csv(metadata,    file.path(out_dir, 'sim_metadata.csv'))
write_csv(truth_dates, file.path(out_dir, 'sim_truth_dates.csv'))
write_csv(truth_pos,   file.path(out_dir, 'sim_truth_positions.csv'))
saveRDS(p, file.path(out_dir, 'sim_settings.rds'))


# =============================================================================
# 4. CHECKS ----
# =============================================================================

# Observation error (km), outside and inside the equinox windows
err_summary <- truth_pos %>%
  left_join(metadata %>% select(geo_trip, error), by = 'geo_trip') %>%
  filter(error != 'none') %>%
  mutate(eq = (abs(as.numeric(difftime(Date_Time, as.POSIXct(paste0(format(Date_Time, '%Y'), '-09-22'), tz = 'UTC'), units = 'days'))) <= p$error$gls$eq_days) |
              (abs(as.numeric(difftime(Date_Time, as.POSIXct(paste0(format(Date_Time, '%Y'), '-03-20'), tz = 'UTC'), units = 'days'))) <= p$error$gls$eq_days)) %>%
  group_by(error, equinox_window = eq) %>%
  summarise(n = n(), median_km = round(median(error_km)), mean_km = round(mean(error_km)),
            p90_km = round(quantile(error_km, 0.9)), .groups = 'drop')
print(err_summary)
write_csv(err_summary, file.path(out_dir, 'checks', 'error_summary.csv'))

# Maps: one example track per movement and error, for each distance class
world <- ne_countries(scale = 'medium', returnclass = 'sf')
for (d in c('short', 'medium', 'long')) {
  ex <- scenarios %>% filter(distance == d, rep == 1) %>% pull(sim_id)
  dat <- positions %>% filter(geo_trip %in% ex) %>%
    left_join(truth_pos %>% select(geo_trip, Date_Time, stage, lon_true, lat_true),
              by = c('geo_trip', 'Date_Time')) %>%
    left_join(metadata %>% select(geo_trip, movement, error), by = 'geo_trip')
  g <- ggplot() +
    geom_sf(data = world, fill = 'grey85', colour = NA) +
    geom_path(data = dat, aes(lon_true, lat_true), colour = 'grey40', linewidth = 0.3) +
    geom_point(data = dat, aes(lon, lat, colour = stage), size = 0.6) +
    geom_point(data = tibble(lon = p$col_lon, lat = p$col_lat), aes(lon, lat),
               shape = 18, size = 3) +
    scale_colour_manual(values = c(B = '#2a6f97', mig = '#e09f3e', W = '#9e2a2b')) +
    coord_sf(xlim = range(dat$lon) + c(-5, 5), ylim = range(dat$lat) + c(-5, 5)) +
    facet_grid(error ~ movement) +
    labs(title = paste('Example tracks:', d, 'distance'),
         subtitle = 'Grey line: true track; points: observed positions (colour = true class)',
         x = NULL, y = NULL, colour = 'True class') +
    theme_bw(base_size = 9)
  ggsave(file.path(out_dir, 'checks', paste0('examples_', d, '.png')), g,
         width = 26, height = 18, units = 'cm', dpi = 200)
}
cat('Done:', nrow(metadata), 'tracks,', nrow(positions), 'positions\n')
