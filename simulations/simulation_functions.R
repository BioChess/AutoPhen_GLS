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

  # --- Dates (true phenology) -----------------------------------------------------
  # Every track has the two equinoxes (+-eq_days around 22 September and 20
  # March). The scenarios only differ in WHERE the equinoxes fall, with the
  # same duration ranges of all stages:
  #   'none', 'gls': both equinoxes during the wintering period (no
  #             migration overlaps an equinox)
  #   'gls_eq': one equinox (chosen at random) at a migration: the outward
  #             migration is centred on the September equinox, or the return
  #             migration is centred on the March equinox (the rest of that
  #             equinox window falls in the stages before and after the
  #             migration); the other equinox is during the wintering period
  # Distances, routes, speeds and durations are drawn again until the dates
  # fit.
  tm <- p$timing
  eq_overlap <- if (sc$error == 'gls_eq') sample(c('september', 'march'), 1) else 'none'
  start <- as.POSIXct(p$start_date, tz = 'UTC')
  yr1 <- as.numeric(format(start, '%Y'))
  rdate <- function(r, y) as.POSIXct(paste0(y, '-', r[1]), tz = 'UTC') +
    runif(1, 0, as.numeric(difftime(as.POSIXct(paste0(y, '-', r[2]), tz = 'UTC'),
                                    as.POSIXct(paste0(y, '-', r[1]), tz = 'UTC'),
                                    units = 'secs')))
  for (k in 1:200) {
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

    # Loop migration (as in the original simulations): the outward and return
    # routes pass through two different waypoints. One waypoint is at a distance
    # D/k (k = 2-4) from the breeding site and the other at D/k from the
    # wintering area, 36-60 degrees to one side of the straight line. If both
    # are on the same side the track has an 8-shape; if not, a round shape.
    side <- function() sample(c(-1, 1), 1) * runif(1, p$waypoint_angle_deg[1],
                                                  p$waypoint_angle_deg[2])
    P_B <- move_point(0, 0, dist_w / sample(p$waypoint_partition, 1), bearing_w + side())
    P_W <- move_point(W2[['x']], W2[['y']], dist_w / sample(p$waypoint_partition, 1),
                      bearing_w + 180 + side())
    if (runif(1) < 0.5) { P_out <- P_B; P_in <- P_W } else { P_out <- P_W; P_in <- P_B }
    if (sc$movement == 'two_W') P_out <- P_B   # P_W is near the second area
    seg_len <- function(a, b) sqrt(sum((b - a)^2))
    O <- c(x = 0, y = 0)

    # outward route: breeding site -> waypoint (stopover) -> wintering area
    n_out1 <- max(1, day(seg_len(O, P_out) / speed_out))
    n_out2 <- max(1, day(seg_len(P_out, W1) / speed_out))
    n_out  <- n_out1 + n_out2
    n_stop <- if (sc$movement == 'stopover') day(runif(1, p$stopover_days[1],
                                                      p$stopover_days[2])) else 0
    # return route: (last) wintering area -> waypoint -> breeding site
    n_in1 <- max(1, day(seg_len(W2, P_in) / speed_in))
    n_in2 <- max(1, day(seg_len(P_in, O) / speed_in))
    n_in  <- n_in1 + n_in2
    n_w2move <- if (sc$movement == 'two_W') {
      day(sqrt(sum((W2 - W1)^2)) / speed_out)
    } else 0
    n_w   <- day(runif(1, p$wintering_days[1], p$wintering_days[2]))
    n_out_all <- n_out + n_stop                 # outward migration (+ stopover)
    n_w_all   <- n_w + n_w2move                 # wintering period (+ move W1-W2)
    step_of <- function(d) round(as.numeric(difftime(d, start, units = 'days')) * ppd)
    eq_win <- function(md, y) {
      e <- step_of(as.POSIXct(paste0(y, '-', md), tz = 'UTC'))
      c(e - day(p$error$gls$eq_days), e + day(p$error$gls$eq_days))
    }
    sep_win <- eq_win('09-22', yr1); mar_win <- eq_win('03-20', yr1 + 1)
    overlaps <- function(a, b, w) a <= w[2] && b >= w[1]
    if (eq_overlap == 'march') {
      # return migration centred on the March equinox
      mid    <- step_of(rdate(tm$march_centre, yr1 + 1))
      t_wdep <- mid - round(n_in / 2)
      t_warr <- t_wdep - n_w_all
      t_dep  <- t_warr - n_out_all
    } else if (eq_overlap == 'september') {
      # outward migration centred on the September equinox
      mid    <- step_of(rdate(tm$september_centre, yr1))
      t_dep  <- mid - round(n_out_all / 2)
      t_warr <- t_dep + n_out_all
      t_wdep <- t_warr + n_w_all
    } else {
      # arrival at the wintering area before the September equinox window
      t_warr <- step_of(rdate(tm$wintering_arrival, yr1))
      t_dep  <- t_warr - n_out_all
      t_wdep <- t_warr + n_w_all
    }
    t_arr <- t_wdep + n_in
    out_eq  <- overlaps(t_dep, t_warr, sep_win)
    back_eq <- overlaps(t_wdep, t_arr, mar_win)
    ok <- t_dep >= day(tm$min_breeding_days) &&
      switch(eq_overlap,
             none      = t_warr < sep_win[1] && t_wdep > mar_win[2],
             september = out_eq && t_wdep > mar_win[2],
             march     = t_warr < sep_win[1] && back_eq)
    if (ok) break
  }
  if (!ok) stop('The migrations do not fit the timing of scenario ', sc$sim_id)
  n_end <- day(runif(1, p$after_arrival_days[1], p$after_arrival_days[2]))
  trips <- sc$movement == 'foraging'

  # --- Legs of the centre of attraction -------------------------------------------
  legs <- breeding_legs(t_dep, 'B1', trips, p)
  # outward migration; the stopover (if any) is at the waypoint. A stopover
  # is a short stationary stage outside the breeding site: for the workflow
  # it is a (short) wintering stage, W1, and the main wintering area is W2.
  # Calling it stopover or wintering only depends on its duration.
  legs <- c(legs, list(
    list(type = 'move', x = P_out[['x']], y = P_out[['y']], n = n_out1, stage = 'mig', label = 'mig1')),
    if (n_stop > 0) list(list(type = 'stay', x = P_out[['x']], y = P_out[['y']],
                              n = n_stop, stage = 'W', label = 'stopover')),
    list(list(type = 'move', x = W1[['x']], y = W1[['y']], n = n_out2, stage = 'mig',
              label = 'mig1')))
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
  legs <- c(legs, list(
    list(type = 'move', x = P_in[['x']], y = P_in[['y']], n = n_in1, stage = 'mig', label = 'mig2'),
    list(type = 'move', x = 0, y = 0, n = n_in2, stage = 'mig', label = 'mig2')),
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
    error = sc$error, seed = sc$seed, equinox_overlap = eq_overlap,
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


# =============================================================================
# 5. LABELS AND PLOTS ----
# =============================================================================
# Same colours and point style as the maps of the workflow
# (functions/segmentation_classification_functions.R)

class_labels    <- c(B = 'Breeding', mig = 'Migration', W = 'Wintering')
class_colours   <- c(Breeding = 'orange2', Migration = 'grey55', Wintering = 'steelblue')
movement_labels <- c(base = 'Base', foraging = 'Foraging trips',
                     stopover = 'Stopover (5-15 days)', two_W = 'Two wintering areas')
error_labels    <- c(none = 'No error (GPS-like)', gls = 'GLS-like error',
                     gls_eq = 'GLS-like error,\nmigration at equinoxes')
distance_labels <- c(short = 'Short (800-2000 km)', medium = 'Medium (3000-5000 km)',
                     long = 'Long (6000-9000 km)')

# Readable labels for the scenario columns
label_scenarios <- function(d) {
  d %>% dplyr::mutate(
    dplyr::across(dplyr::any_of('movement'),
                  ~ factor(movement_labels[.x], levels = movement_labels)),
    dplyr::across(dplyr::any_of('error'),
                  ~ factor(error_labels[.x], levels = error_labels)),
    dplyr::across(dplyr::any_of('distance'),
                  ~ factor(distance_labels[.x], levels = distance_labels)))
}

theme_sim <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_line(colour = 'grey92'),
                 panel.border = ggplot2::element_rect(colour = 'black', fill = NA),
                 strip.background = ggplot2::element_rect(fill = 'white', colour = 'black'),
                 strip.text = ggplot2::element_text(face = 'bold'),
                 legend.position = 'bottom',
                 legend.title = ggplot2::element_text(face = 'bold'))

# Example tracks: rows = observation error, columns = movement type.
# dat: rows of sim_truth_positions.csv (+ movement and error), one track per
# panel.
plot_sim_tracks <- function(dat, title = NULL) {
  dat <- dat %>%
    label_scenarios() %>%
    dplyr::mutate(Class = factor(class_labels[stage], levels = class_labels))
  ggplot2::ggplot(dat) +
    ggplot2::geom_path(ggplot2::aes(x_true_km, y_true_km), colour = 'grey20',
                       linewidth = 0.25, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x_km, y_km, fill = Class), shape = 21,
                        colour = 'black', stroke = 0.15, size = 1.3) +
    ggplot2::annotate('point', x = 0, y = 0, shape = 23, size = 2.6,
                      fill = NA, colour = 'black', stroke = 0.9) +
    ggplot2::scale_fill_manual(values = class_colours, name = 'True class',
                               drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ',')) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ',')) +
    ggplot2::facet_grid(error ~ movement) +
    ggplot2::labs(title = title, x = 'x (km)', y = 'y (km)',
                  caption = paste('Line: true track. Points: observed positions,',
                                  'coloured by true class. Diamond: breeding site.')) +
    theme_sim
}
