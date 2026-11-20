# =============================================================================
# SIMULATIONS 01: SIMULATED MIGRATORY TRACKS AT REAL SCALES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Simulates tracks of a central-place forager that breeds at a
#              breeding site, migrates to one or two wintering areas and comes
#              back. The tracks are not linked to any place: they are simulated
#              on a plane in km, with the breeding site at the origin (0, 0)
#              and the wintering area in a random direction. Movement is
#              simulated with marcher (Gurarie et al. 2017) around a moving
#              centre of attraction, with time in days. Observation error of
#              different tracking methods is added to the true positions.
#
#   Scenarios (full factorial, n_rep tracks each):
#     distance  short (800-2000 km), medium (3000-5000 km), long (6000-9000 km)
#     movement  base      one breeding area and one wintering area
#               foraging  central-place foraging trips (4-20 days) from the
#                         breeding site during the breeding periods
#               stopover  a stopover during the outward migration
#               two_W     two wintering areas
#     error     none      no observation error (e.g. GPS)
#               gls       GLS-like error (larger in the north-south axis,
#                         autocorrelated, north-south error inflated around
#                         the equinoxes, missing positions); migrations away
#                         from the equinoxes
#               gls_eq    the same error, with migrations overlapping the
#                         equinoxes
#
#   Coordinates for the workflow: the workflow needs longitude and latitude.
#   The plane (km) is placed on an arbitrary reference point (ref_lon,
#   ref_lat) with an azimuthal equidistant projection. The workflow projects
#   the tracks again with the same projection centred on the breeding site, so
#   it works with the same km coordinates and the results do not depend on
#   the reference point. Only the dates have a real meaning (equinoxes).
#
# Output (simulations/data/):
#   sim_positions.csv        observed positions (geo_trip, Date_Time, lon, lat):
#                            the same format as data/All_positions.csv
#   sim_metadata.csv         one row per track (format of data/All_metadata.csv
#                            + scenario)
#   sim_truth_dates.csv      true dates and parameters of each track
#   sim_truth_positions.csv  observed and true positions in km, and true class
#                            (B, mig, W) at each observed time
#   checks/                  plots of example tracks and summary of the error
#
# Run from the project root. Needs the package marcher (CRAN).
# =============================================================================


# =============================================================================
# 1. PACKAGES AND SETTINGS ----
# =============================================================================

pacman::p_load('tidyverse', 'sf', 'marcher')
source('simulations/simulation_functions.R')

out_dir <- file.path('simulations', 'data')
dir.create(file.path(out_dir, 'checks'), recursive = TRUE, showWarnings = FALSE)

p <- list(
  positions_per_day = 2,
  start_date = '2021-07-01',          # first position (breeding period)
  # Reference point to convert km to longitude/latitude (arbitrary; see above)
  ref_lon = 0, ref_lat = 0,

  # Number of tracks per scenario
  n_rep = 50,

  # --- Movement --------------------------------------------------------------------
  # Distance from the breeding site to the (first) wintering area, by class (km)
  distance_km = list(short = c(800, 2000), medium = c(3000, 5000),
                     long = c(6000, 9000)),
  bearing_deg = c(0, 360),            # direction of the wintering area (random)
  speed_km_day = c(300, 600),         # migration speed
  # loop migration: waypoints of the outward and return routes (see
  # simulation_functions.R)
  waypoint_partition = 2:4,           # waypoint at D/k from one end (k = 2-4)
  waypoint_angle_deg = c(36, 60),     # angle from the straight line (pi/5-pi/3)
  range_radius_km = c(100, 250),      # radius with 95% of the positions around
                                      # the centre (marcher A = pi * r^2)
  tau_z_days = c(2.5, 3.5),           # position autocorrelation (marcher tau.z)
  tau_v_days = c(0.02, 0.03),         # velocity autocorrelation (marcher tau.v)
  min_wintering_days = 30,            # minimum stay in each wintering area
  after_arrival_days = c(15, 45),     # breeding period after the return
  # foraging trips (movement 'foraging')
  trip_gap_days = c(1, 15),           # time at the breeding site between trips
  trip_days = c(4, 20),               # duration of a trip
  trip_km = c(150, 1000),             # maximum distance from the breeding site
  # stopover (movement 'stopover'): at the waypoint of the outward route
  stopover_days = c(5, 15),
  # second wintering area (movement 'two_W')
  w2_km = c(500, 1500),               # distance between the two wintering areas

  # --- Timing (month-day ranges) ----------------------------------------------------
  # normal:  migrations away from the equinoxes (+-20 days)
  # equinox: outward and return migrations overlap the equinoxes
  timing = list(normal  = list(B1_dep = c('10-20', '11-15'), B2_arr = c('02-05', '02-25')),
                equinox = list(B1_dep = c('09-10', '10-02'), B2_arr = c('03-08', '03-28'))),

  # --- Observation error ----------------------------------------------------------
  error = list(
    # GLS: mean error of about 150-200 km, larger in latitude (north-south
    # axis) (Phillips et al. 2004; Halpin et al. 2021); errors correlated
    # between consecutive positions; north-south position almost
    # uninformative within +-eq_days of each equinox
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

  # The first and last two positions are set at the breeding site (known
  # deployment and recovery), as in the real data
  n <- nrow(obs)
  obs$x[c(1:2, (n - 1):n)] <- 0; obs$y[c(1:2, (n - 1):n)] <- 0

  ll <- km_to_lonlat(obs$x, obs$y, p$ref_lon, p$ref_lat)
  positions[[i]] <- tibble(geo_trip = sc$sim_id, Date_Time = obs$Date_Time,
                           lon = ll$lon, lat = ll$lat)
  truth_pos[[i]] <- tibble(geo_trip = sc$sim_id, Date_Time = obs$Date_Time,
                           x_km = obs$x, y_km = obs$y,
                           x_true_km = obs$x_true, y_true_km = obs$y_true,
                           error_km = sqrt((obs$x - obs$x_true)^2 + (obs$y - obs$y_true)^2),
                           stage = obs$stage, leg = obs$leg)
  truth_dates[[i]] <- s$dates
  if (i %% 50 == 0) cat('  ', i, 'tracks\n')
}
positions   <- bind_rows(positions)
truth_pos   <- bind_rows(truth_pos)
truth_dates <- bind_rows(truth_dates)

# The breeding site is the reference point (0, 0 km)
metadata <- positions %>%
  group_by(geo_trip) %>%
  summarise(n_positions = n(), start_date = min(Date_Time), end_date = max(Date_Time),
            .groups = 'drop') %>%
  mutate(Ring = geo_trip, Sp = 'SIM', Species = 'Simulated',
         Common_name = 'Simulated track', Colony = 'Breeding site',
         Colony.lat = p$ref_lat, Colony.lon = p$ref_lon,
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
eq_dist <- function(d, md) {
  abs(as.numeric(difftime(d, as.POSIXct(paste0(format(d, '%Y'), md), tz = 'UTC'),
                          units = 'days')))
}
err_summary <- truth_pos %>%
  left_join(metadata %>% select(geo_trip, error), by = 'geo_trip') %>%
  filter(error != 'none') %>%
  mutate(eq = eq_dist(Date_Time, '-09-22') <= p$error$gls$eq_days |
              eq_dist(Date_Time, '-03-20') <= p$error$gls$eq_days) %>%
  group_by(error, equinox_window = eq) %>%
  summarise(n = n(), median_km = round(median(error_km)), mean_km = round(mean(error_km)),
            p90_km = round(quantile(error_km, 0.9)), .groups = 'drop')
print(err_summary)
write_csv(err_summary, file.path(out_dir, 'checks', 'error_summary.csv'))

# Plots: one example track per movement and error, for each distance class
# (km from the breeding site; no map, the tracks are not linked to any place)
for (d in c('short', 'medium', 'long')) {
  ex <- scenarios %>% filter(distance == d, rep == 1) %>% pull(sim_id)
  dat <- truth_pos %>% filter(geo_trip %in% ex) %>%
    left_join(metadata %>% select(geo_trip, movement, error), by = 'geo_trip')
  g <- plot_sim_tracks(dat, title = paste('Example tracks:', distance_labels[[d]]))
  ggsave(file.path(out_dir, 'checks', paste0('examples_', d, '.png')), g,
         width = 26, height = 19, units = 'cm', dpi = 200)
}
cat('Done:', nrow(metadata), 'tracks,', nrow(positions), 'positions\n')
