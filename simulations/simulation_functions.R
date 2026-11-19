# =============================================================================
# SIMULATIONS: FUNCTIONS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Functions to simulate migratory tracks at real scales (km and
#              days) with the package marcher, and to add the observation error
#              of different tracking methods. Used by 01_simulate_tracks.R.
#
# How a track is built:
#   1. The ATTRACTION CENTRE (mu) of the animal moves between stationary
#      areas: breeding site -> (stopover) -> wintering area(s) -> breeding
#      site. During the breeding periods the centre can also make
#      out-and-back excursions from the breeding site (central-place foraging
#      trips).
#   2. marcher::simulate_shift() simulates the movement of the animal around
#      that centre (Ornstein-Uhlenbeck process with position and velocity
#      autocorrelation; Gurarie et al. 2017).
#   3. Observation error is added to the true positions (none, PTT-like or
#      GLS-like, with or without equinox effects and data gaps).
#   4. Positions are on a plane in km, with the breeding site at (0, 0). For
#      the workflow, they are converted to longitude and latitude with an
#      azimuthal equidistant projection centred on an arbitrary reference
#      point (see 01_simulate_tracks.R).
# =============================================================================


# =============================================================================
# 1. CENTRE OF ATTRACTION (mu) ----
# =============================================================================

# Point at a distance (km) and bearing (degrees from north) from (x, y)
move_point <- function(x, y, dist, bearing) {
  b <- bearing * pi / 180
  c(x = x + dist * sin(b), y = y + dist * cos(b))
}

# Builds the path of the centre of attraction from a list of 'legs'. Each leg
# is the centre staying at a point (type 'stay') or moving in a straight line
# to a point (type 'move'), during n steps. 'stage' is the true class of the
# positions of that leg (B, mig or W).
legs_to_mu <- function(legs) {
  out <- vector('list', length(legs))
  x0 <- legs[[1]]$x; y0 <- legs[[1]]$y
  for (i in seq_along(legs)) {
    l <- legs[[i]]
    if (l$type == 'stay') {
      xs <- rep(l$x, l$n); ys <- rep(l$y, l$n)
    } else {
      f <- seq_len(l$n) / l$n
      xs <- x0 + (l$x - x0) * f; ys <- y0 + (l$y - y0) * f
    }
    out[[i]] <- data.frame(mu_x = xs, mu_y = ys, stage = l$stage,
                           leg = l$label)
    x0 <- l$x; y0 <- l$y
  }
  do.call(rbind, out)
}

# Breeding period at the breeding site, with or without foraging trips
breeding_legs <- function(n, label, trips, p) {
  if (!trips) return(list(list(type = 'stay', x = 0, y = 0, n = n,
                               stage = 'B', label = label)))
  legs <- list(); used <- 0
  while (used < n) {
    # time at the breeding site before the next trip
    stay <- min(n - used, round(runif(1, p$trip_gap_days[1], p$trip_gap_days[2]) *
                                  p$positions_per_day))
    legs[[length(legs) + 1]] <- list(type = 'stay', x = 0, y = 0, n = max(1, stay),
                                     stage = 'B', label = label)
    used <- used + max(1, stay)
    # trip: the centre goes out and comes back (half of the duration each way)
    dur <- round(runif(1, p$trip_days[1], p$trip_days[2]) * p$positions_per_day)
    if (used + dur > n - 2) break
    far <- move_point(0, 0, runif(1, p$trip_km[1], p$trip_km[2]), runif(1, 0, 360))
    half <- max(1, round(dur / 2))
    legs[[length(legs) + 1]] <- list(type = 'move', x = far[['x']], y = far[['y']],
                                     n = half, stage = 'B', label = label)
    legs[[length(legs) + 1]] <- list(type = 'move', x = 0, y = 0, n = dur - half,
                                     stage = 'B', label = label)
    used <- used + dur
  }
  if (used < n) legs[[length(legs) + 1]] <- list(type = 'stay', x = 0, y = 0,
                                                 n = n - used, stage = 'B',
                                                 label = label)
  legs
}


# =============================================================================
# 2. ONE SIMULATED TRACK ----
# =============================================================================
# sc: one row of the scenario table (distance class, movement, error, seed)
# p:  list of settings (01_simulate_tracks.R)
# Returns the true positions (km) with their true class, and the true dates.

simulate_track <- function(sc, p) {
  set.seed(sc$seed)
  ppd <- p$positions_per_day
  day <- function(d) round(d * ppd)               # days -> steps

  # --- Distances and directions -------------------------------------------------
  dist_w <- runif(1, p$distance_km[[sc$distance]][1], p$distance_km[[sc$distance]][2])
  bearing_w <- runif(1, p$bearing_deg[1], p$bearing_deg[2])
  W1 <- move_point(0, 0, dist_w, bearing_w)
  if (sc$movement == 'two_W') {
    W2 <- move_point(W1[['x']], W1[['y']], runif(1, p$w2_km[1], p$w2_km[2]),
                     runif(1, 0, 360))
  } else W2 <- W1
  speed_out <- runif(1, p$speed_km_day[1], p$speed_km_day[2])
  speed_in  <- runif(1, p$speed_km_day[1], p$speed_km_day[2])

  # --- Dates (true phenology) -----------------------------------------------------
  # Migration away from the equinoxes ('none', 'gls') or overlapping them
  # ('gls_eq'). Dates are drawn until the wintering stay is long enough.
  timing <- if (sc$error == 'gls_eq') p$timing$equinox else p$timing$normal
  start <- as.POSIXct(p$start_date, tz = 'UTC')
  yr1 <- as.numeric(format(start, '%Y'))
  rdate <- function(r, y) as.POSIXct(paste0(y, '-', r[1]), tz = 'UTC') +
    runif(1, 0, as.numeric(difftime(as.POSIXct(paste0(y, '-', r[2]), tz = 'UTC'),
                                    as.POSIXct(paste0(y, '-', r[1]), tz = 'UTC'),
                                    units = 'secs')))
  for (k in 1:50) {
    dep <- rdate(timing$B1_dep, yr1)
    arr <- rdate(timing$B2_arr, yr1 + 1)
    # outward migration: straight line (+ stopover); return: straight line
    route_out <- dist_w
    n_out <- day(route_out / speed_out)
    n_stop <- if (sc$movement == 'stopover') day(runif(1, p$stopover_days[1],
                                                      p$stopover_days[2])) else 0
    dist_back <- sqrt(W2[['x']]^2 + W2[['y']]^2)
    n_in <- day(dist_back / speed_in)
    n_w2move <- if (sc$movement == 'two_W') {
      day(sqrt(sum((W2 - W1)^2)) / speed_out)
    } else 0
    t_dep <- round(as.numeric(difftime(dep, start, units = 'days')) * ppd)
    t_arr <- round(as.numeric(difftime(arr, start, units = 'days')) * ppd)
    n_w <- t_arr - n_in - (t_dep + n_out + n_stop) - n_w2move
    min_w <- day(p$min_wintering_days) * ifelse(sc$movement == 'two_W', 2, 1)
    if (n_w >= min_w) break
  }
  if (n_w < min_w) stop('Wintering stay too short for scenario ', sc$sim_id)
  n_end <- day(runif(1, p$after_arrival_days[1], p$after_arrival_days[2]))
  trips <- sc$movement == 'foraging'

  # --- Legs of the centre of attraction -------------------------------------------
  legs <- breeding_legs(t_dep, 'B1', trips, p)
  if (sc$movement == 'stopover') {
    f <- runif(1, p$stopover_position[1], p$stopover_position[2])
    S <- move_point(0, 0, dist_w * f, bearing_w + runif(1, -20, 20))
    n1 <- max(1, round(n_out * f))
    legs <- c(legs, list(
      list(type = 'move', x = S[['x']], y = S[['y']], n = n1, stage = 'mig', label = 'mig1'),
      list(type = 'stay', x = S[['x']], y = S[['y']], n = n_stop, stage = 'mig', label = 'stopover'),
      list(type = 'move', x = W1[['x']], y = W1[['y']], n = max(1, n_out - n1),
           stage = 'mig', label = 'mig1')))
  } else {
    legs <- c(legs, list(list(type = 'move', x = W1[['x']], y = W1[['y']], n = n_out,
                              stage = 'mig', label = 'mig1')))
  }
  if (sc$movement == 'two_W') {
    n_w1 <- round(n_w * runif(1, 0.35, 0.65))
    legs <- c(legs, list(
      list(type = 'stay', x = W1[['x']], y = W1[['y']], n = n_w1, stage = 'W', label = 'W1'),
      list(type = 'move', x = W2[['x']], y = W2[['y']], n = n_w2move, stage = 'mig', label = 'mig_W'),
      list(type = 'stay', x = W2[['x']], y = W2[['y']], n = n_w - n_w1, stage = 'W', label = 'W2')))
  } else {
    legs <- c(legs, list(list(type = 'stay', x = W1[['x']], y = W1[['y']], n = n_w,
                              stage = 'W', label = 'W1')))
  }
  legs <- c(legs, list(list(type = 'move', x = 0, y = 0, n = n_in, stage = 'mig',
                            label = 'mig2')),
            breeding_legs(n_end, 'B2', trips, p))
  mu <- legs_to_mu(legs)

  # --- Movement around the centre (marcher) ------------------------------------------
  # A = area (km2) where 95% of the positions are expected around the centre
  A <- pi * runif(1, p$range_radius_km[1], p$range_radius_km[2])^2
  tau <- c(tau.z = runif(1, p$tau_z_days[1], p$tau_z_days[2]) * ppd,
           tau.v = runif(1, p$tau_v_days[1], p$tau_v_days[2]) * ppd)
  T <- seq_len(nrow(mu))
  sim <- marcher::simulate_shift(T, tau = tau, mu = as.matrix(
    data.frame(x = mu$mu_x, y = mu$mu_y)), A = A)

  true <- data.frame(sim_id = sc$sim_id, step = T,
                     Date_Time = start + (T - 1) * 86400 / ppd,
                     x_true = sim$X, y_true = sim$Y,
                     mu_x = mu$mu_x, mu_y = mu$mu_y,
                     stage = mu$stage, leg = mu$leg)

  # --- True dates ------------------------------------------------------------------
  # Departure = first step of the centre leaving the area; arrival = first step
  # with the centre in the area (the same definition as the KDE dates: first
  # position outside / inside)
  st <- true$stage
  first_mig <- min(which(st == 'mig'))
  first_W   <- min(which(st == 'W'))
  last_W    <- max(which(st == 'W'))
  first_B2  <- min(which(true$leg == 'B2'))
  dates <- data.frame(
    sim_id = sc$sim_id, distance = sc$distance, movement = sc$movement,
    error = sc$error, seed = sc$seed,
    dist_km = round(dist_w), bearing = round(bearing_w),
    speed_out = round(speed_out), speed_in = round(speed_in),
    range_radius_km = round(sqrt(A / pi)), tau_z_days = round(tau[['tau.z']] / ppd, 2),
    tau_v_days = round(tau[['tau.v']] / ppd, 3),
    B1_dep = as.Date(true$Date_Time[first_mig]),
    W_arr  = as.Date(true$Date_Time[first_W]),
    W_dep  = as.Date(true$Date_Time[last_W + 1]),
    B2_arr = as.Date(true$Date_Time[first_B2]))
  list(true = true, dates = dates)
}


# =============================================================================
# 3. OBSERVATION ERROR ----
# =============================================================================

# AR(1) series with standard deviation sd and correlation phi between
# consecutive positions
ar1 <- function(n, sd, phi) {
  e <- numeric(n); e[1] <- rnorm(1, 0, sd)
  for (i in 2:n) e[i] <- phi * e[i - 1] + rnorm(1, 0, sd * sqrt(1 - phi^2))
  e
}

# Adds the error of one tracking method to the true positions (km)
#   'none' : no error (e.g. GPS)
#   'gls', 'gls_eq' : GLS-like error, larger in the north-south axis (y,
#            latitude in GLS), autocorrelated; around each equinox the y error
#            is strongly inflated; random missing positions and short gaps
#   'ptt'  : Argos-like error (optional block)
add_error <- function(true, error, p, seed) {
  set.seed(seed + 100000)
  n <- nrow(true)
  obs <- true; obs$x <- true$x_true; obs$y <- true$y_true
  if (error == 'none') return(obs)

  if (error == 'ptt') {
    e <- p$error$ptt
    obs$x <- obs$x + rnorm(n, 0, e$sd_km); obs$y <- obs$y + rnorm(n, 0, e$sd_km)
    return(obs)
  }

  e <- p$error$gls
  ex <- ar1(n, e$sd_lon_km, e$phi)
  ey <- ar1(n, e$sd_lat_km, e$phi)
  # Equinoxes: y (latitude) error inflated within +-eq_days of each equinox
  doy_dist <- function(d, m, dd) {
    eq <- as.POSIXct(paste0(format(d, '%Y'), '-', m, '-', dd), tz = 'UTC')
    abs(as.numeric(difftime(d, eq, units = 'days')))
  }
  near_eq <- doy_dist(true$Date_Time, '09', '22') <= e$eq_days |
             doy_dist(true$Date_Time, '03', '20') <= e$eq_days
  ey_eq <- ar1(n, e$sd_lat_km * e$eq_lat_factor, e$eq_phi)
  ey[near_eq] <- ey_eq[near_eq]
  obs$x <- obs$x + ex; obs$y <- obs$y + ey
  obs$equinox_true <- near_eq

  # Missing data: random positions and short gaps (not in the first and last
  # positions, which are fixed at the breeding site)
  keep <- runif(n) > e$p_missing
  n_gaps <- rpois(1, e$gaps_per_track)
  for (g in seq_len(n_gaps)) {
    s <- sample(5:(n - 10), 1)
    keep[s:min(n - 5, s + round(runif(1, e$gap_days[1], e$gap_days[2]) *
                                   p$positions_per_day))] <- FALSE
  }
  keep[c(1:2, (n - 1):n)] <- TRUE
  obs[keep, ]
}


# =============================================================================
# 4. KM TO LONGITUDE / LATITUDE ----
# =============================================================================

km_to_lonlat <- function(x, y, ref_lon, ref_lat) {
  aeqd <- paste0('+proj=aeqd +lat_0=', ref_lat, ' +lon_0=', ref_lon,
                 ' +x_0=0 +y_0=0 +units=km +datum=WGS84')
  pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c('x', 'y'), crs = aeqd)
  ll <- sf::st_coordinates(sf::st_transform(pts, 4326))
  data.frame(lon = ll[, 1], lat = ll[, 2])
}
