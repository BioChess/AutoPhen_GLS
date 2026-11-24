# =============================================================================
# BASELINE METHODS: FUNCTIONS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Three simple alternatives to the workflow, used to compare it
#              with existing approaches (reviewer comment R2-M7). All methods
#              use the same positions as the workflow (projection centred on
#              the colony and the same cleaning), and the same final rules to
#              build the stages, so only the classification step changes.
#
#   segclust_simple  The same segclust2d segmentation as the workflow, with a
#                    simple rule for each segment: breeding if its mean NSCD
#                    is below the NSCD limit of the workflow (R2n_limit);
#                    otherwise migration if its net speed (straight-line
#                    distance between its first and last positions / duration)
#                    is higher than 'simple_speed_km_day', and wintering if
#                    not. No double clustering, no agreement score and no
#                    correction rules.
#   migrateR         Net squared displacement models (Bunnefeld et al. 2011;
#                    Spitz et al. 2017) fitted to the NSCD (distance to the
#                    colony). If the best model (AIC) is 'migrant' or
#                    'mixmig', the four dates are the start and end of the
#                    outward and return movements (5% and 95% of the change in
#                    NSCD; migrateR::mvmt2dt). Other models: no migration.
#   hmm              Hidden Markov model with three states (momentuHMM;
#                    McClintock & Michelot 2018) on step length and distance to
#                    the colony (both gamma distributions). States are named
#                    by their means: breeding = shortest distance to the
#                    colony; of the other two, migration = longest steps,
#                    wintering = the other one. Classes from the Viterbi
#                    sequence. The model is fitted from several starting
#                    values and the best fit is kept.
#
#   Common final rules (as in the workflow, step 03):
#     1. Short stages (fewer than lmin positions) between two stages of the
#        same class take that class (check_phen()).
#     2. Movements between two wintering stages, without breeding between
#        them, are part of the wintering period (W-mig-W = W).
#   Dates from the classes (first day of each change):
#     B1_dep = first position after the first breeding stage
#     W_arr  = first wintering position
#     W_dep  = first position after the last wintering position
#     B2_arr = first breeding position after the last wintering position
#   migrateR gives its own dates (above); its classes are built from them.
# =============================================================================


# =============================================================================
# 1. POSITIONS (as in classify_track()) ----
# =============================================================================

# Projection centred on the colony (km), trajectory and the same cleaning as
# the workflow. Returns the cleaned positions (date, x, y in km, R2n in km^2).
prepare_track <- function(track, info, p, env) {
  track <- dplyr::arrange(track, Date_Time)
  proj_colony <- paste0('+proj=aeqd +lat_0=', info$Colony.lat,
                        ' +lon_0=', info$Colony.lon, ' +units=km +datum=WGS84')
  xy_km <- sf::sf_project(from = sf::st_crs(4326), to = proj_colony,
                          pts = as.matrix(track[, c('lon', 'lat')]))
  ltraj <- adehabitatLT::as.ltraj(xy = xy_km, date = track$Date_Time,
                                  id = info$geo_trip)
  df <- ltraj[[1]]
  pos <- match(df$date, track$Date_Time)
  df$lon <- track$lon[pos]; df$lat <- track$lat[pos]
  df$R2n <- df$x^2 + df$y^2

  next_date <- function(start, month_day) {
    d <- as.Date(paste0(lubridate::year(start), month_day))
    if (d < as.Date(start)) d <- as.Date(paste0(lubridate::year(start) + 1, month_day))
    d
  }
  Sep_eq <- next_date(track$Date_Time[1], '-09-21')
  Mar_eq <- next_date(track$Date_Time[1], '-03-21')
  eq_days <- if (p$use_equinox) p$eq_window_days else NULL

  df2 <- env$clean_positions(data = df, dist_threshold = p$adv$max_step)
  invisible(utils::capture.output(
    df3 <- env$remove_duplicates(data = df2, equinox_days = eq_days,
                                 Sep_eq = Sep_eq, Mar_eq = Mar_eq,
                                 dup_threshold = p$adv$max_repeated,
                                 edge_positions = p$adv$edge_positions)))
  df3
}


# =============================================================================
# 2. COMMON FINAL RULES AND DATES ----
# =============================================================================

# Classes B / mig / W of each position -> final classes (rules 1 and 2)
final_classes <- function(cls, lmin, env) {
  d <- data.frame(Breed2 = cls, stringsAsFactors = FALSE)
  d <- env$check_phen(d, 'Breed2', min_length = lmin)
  r <- rle(d$Breed2)
  is_W <- r$values == 'W'
  for (k in seq_along(r$values)) {
    if (r$values[k] == 'mig') {
      prev_W <- which(is_W[seq_len(k - 1)])
      next_W <- which(is_W[-seq_len(k)]) + k
      if (length(prev_W) > 0 && length(next_W) > 0) {
        between <- r$values[max(prev_W):min(next_W)]
        if (!any(between == 'B')) r$values[k] <- 'W'
      }
    }
  }
  inverse.rle(r)
}

# Four dates from the classes of a track (dates of the positions, as Date)
dates_from_classes <- function(date, cls) {
  na_date <- as.Date(NA)
  out <- list(B1_dep = na_date, W_arr = na_date, W_dep = na_date, B2_arr = na_date)
  n <- length(cls)
  if (n == 0) return(out)
  d <- as.Date(date)
  # First breeding stage (the track starts at the colony)
  first_B <- which(cls == 'B')[1]
  if (!is.na(first_B)) {
    end_B1 <- first_B
    while (end_B1 < n && cls[end_B1 + 1] == 'B') end_B1 <- end_B1 + 1
    if (end_B1 < n) out$B1_dep <- d[end_B1 + 1]
  }
  w <- which(cls == 'W')
  if (length(w) > 0) {
    out$W_arr <- d[w[1]]
    last_W <- max(w)
    if (last_W < n) out$W_dep <- d[last_W + 1]
    b2 <- which(cls == 'B' & seq_len(n) > last_W)
    if (length(b2) > 0) out$B2_arr <- d[b2[1]]
  }
  out
}


# =============================================================================
# 3. METHODS ----
# =============================================================================

# --- 3.1. segclust2d with a simple rule ----------------------------------------
# 'mode_segclust': the segclust2d result of the workflow (from the cache) or
# NULL (then it is calculated with the same settings).
baseline_segclust_simple <- function(df3, p, env, cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    mode_segclust <- readRDS(cache_file)
  } else {
    mode_segclust <- segclust2d::segclust(
      df3, Kmax = p$Kmax, lmin = p$lmin, ncluster = p$ncluster,
      seg.var = c('x', 'y'), order.var = 'x',
      diag.var = c('dx', 'dy', 'R2n', 'x', 'y', 'dist'), subsample = FALSE)
  }
  track <- segclust2d::augment(mode_segclust)
  segs  <- segclust2d::segment(mode_segclust)
  limit <- env$R2n_limit(track, plot = FALSE, Kmax = p$adv$nscd_Kmax,
                         lmin = p$lmin, ncluster = p$adv$nscd_ncluster)
  seg_class <- sapply(seq_len(nrow(segs)), function(i) {
    rows <- segs$begin[i]:segs$end[i]
    if (mean(track$R2n[rows]) < limit) return('B')
    days <- as.numeric(difftime(track$date[max(rows)], track$date[min(rows)],
                                units = 'days'))
    net  <- sqrt((track$x[max(rows)] - track$x[min(rows)])^2 +
                   (track$y[max(rows)] - track$y[min(rows)])^2)
    speed <- net / max(days, 0.5)
    if (speed > p$simple_speed_km_day) 'mig' else 'W'
  })
  cls <- seg_class[findInterval(seq_len(nrow(track)), segs$begin)]
  data.frame(date = track$date, class = cls, stringsAsFactors = FALSE)
}

# --- 3.2. migrateR ---------------------------------------------------------------
# NSD models on the distance to the colony (NSCD). migrateR reads the NSD from
# the 'infolocs' of the trajectory (km^2). The models are fitted from several
# starting values (time of the outward movement and of the stay) and the fit
# with the lowest AIC of each model is kept.
baseline_migrater <- function(df3, p) {
  dd <- as.numeric(difftime(df3$date, df3$date[1], units = 'days'))
  span <- max(dd)
  lt <- adehabitatLT::as.ltraj(xy = df3[, c('x', 'y')] * 1000, date = df3$date,
                               id = 'track',
                               infolocs = data.frame(nsd = df3$R2n))
  fits <- list()
  for (t0 in p$migrater_theta_start * span) {
    for (r0 in p$migrater_rho_start * span) {
      pe <- migrateR::pEst(s.t = t0, u.t = span, s.r = r0, u.r = span,
                           u.p = p$migrater_phi_max, u.p2 = p$migrater_phi_max,
                           s.p = min(7, p$migrater_phi_max / 2),
                           s.p2 = min(7, p$migrater_phi_max / 2))
      m <- tryCatch(suppressWarnings(migrateR::mvmtClass(lt, fam = 'nsd',
                                                          p.est = pe)),
                    error = function(e) NULL)
      if (is.null(m)) next
      for (mod in names(m[[1]]@models)) {
        a <- tryCatch(AIC(m[[1]]@models[[mod]]), error = function(e) NA)
        if (is.na(a)) next
        if (is.null(fits[[mod]]) || a < fits[[mod]]$aic) {
          fits[[mod]] <- list(aic = a, mvmt = m)
        }
      }
    }
  }
  if (length(fits) == 0) stop('migrateR: no model could be fitted')
  aics <- sapply(fits, `[[`, 'aic')
  best <- names(which.min(aics))
  out <- data.frame(date = df3$date, class = 'B', stringsAsFactors = FALSE)
  dates <- list(B1_dep = as.Date(NA), W_arr = as.Date(NA),
                W_dep = as.Date(NA), B2_arr = as.Date(NA))
  if (best %in% c('migrant', 'mixmig')) {
    # Dates as in migrateR::mvmt2dt() (start and end of each movement = when
    # the NSCD has changed by p and 1 - p of its total change). They are
    # calculated here from the coefficients because mvmt2dt() fails in
    # recent R versions (theta2() compares class() with '==').
    cf <- stats::coef(fits[[best]]$mvmt[[1]]@models[[best]])
    lq <- log((1 - p$migrater_p) / p$migrater_p)
    theta2 <- cf[['theta']] + 2 * cf[['phi']] + cf[['rho']] + 2 * cf[['phi2']]
    dday <- c(cf[['theta']] - cf[['phi']] * lq, cf[['theta']] + cf[['phi']] * lq,
              theta2 - cf[['phi2']] * lq, theta2 + cf[['phi2']] * lq)
    # dday = days since the first position (migrateR without 'stdt')
    ev <- df3$date[1] + dday * 86400
    out$class[df3$date >= ev[1] & df3$date < ev[2]] <- 'mig'
    out$class[df3$date >= ev[2] & df3$date < ev[3]] <- 'W'
    out$class[df3$date >= ev[3] & df3$date < ev[4]] <- 'mig'
    dates <- list(B1_dep = as.Date(ev[1]), W_arr = as.Date(ev[2]),
                  W_dep = as.Date(ev[3]), B2_arr = as.Date(ev[4]))
  }
  list(classes = out, dates = dates, model = best,
       aic = paste(names(aics), round(aics, 1), sep = '=', collapse = ';'))
}

# --- 3.3. Hidden Markov model -----------------------------------------------------
baseline_hmm <- function(df3, p) {
  d <- data.frame(ID = 'track', x = df3$x, y = df3$y,
                  dist_col = pmax(sqrt(df3$R2n), 1))
  d <- momentuHMM::prepData(d, type = 'UTM', coordNames = c('x', 'y'),
                            covNames = 'dist_col')
  d$dist_col <- pmax(d$dist_col, 1)
  st <- d$step[!is.na(d$step)]
  zero_step <- any(st == 0)
  q  <- function(v, pr) as.numeric(stats::quantile(v, pr, na.rm = TRUE))
  dc <- d$dist_col

  best <- NULL
  set.seed(1)
  for (k in seq_len(p$hmm_n_starts)) {
    jit <- if (k == 1) 1 else stats::runif(1, 0.6, 1.6)
    # states: 1 = breeding, 2 = migration, 3 = wintering (starting values)
    step_mean <- c(q(st, 0.3), q(st, 0.9), q(st, 0.3)) * jit
    dist_mean <- c(q(dc, 0.1), q(dc, 0.5), q(dc, 0.9)) * c(jit, 1, 1 / jit)
    step_par <- c(step_mean, step_mean * 0.8)
    if (zero_step) step_par <- c(step_par, rep(0.01, 3))
    dist_par <- c(dist_mean, dist_mean * 0.5)
    fit <- tryCatch(suppressMessages(suppressWarnings(
      momentuHMM::fitHMM(d, nbStates = 3,
                         dist = list(step = 'gamma', dist_col = 'gamma'),
                         Par0 = list(step = step_par, dist_col = dist_par),
                         stateNames = c('s1', 's2', 's3')))),
      error = function(e) NULL)
    if (is.null(fit)) next
    if (is.null(best) || fit$mod$minimum < best$mod$minimum) best <- fit
  }
  if (is.null(best)) stop('HMM: the model could not be fitted')

  est <- best$mle
  step_m <- est$step['mean', ]
  dist_m <- est$dist_col['mean', ]
  b_state <- which.min(dist_m)
  others  <- setdiff(1:3, b_state)
  m_state <- others[which.max(step_m[others])]
  w_state <- setdiff(others, m_state)
  vit <- momentuHMM::viterbi(best)
  lab <- character(3); lab[b_state] <- 'B'; lab[m_state] <- 'mig'; lab[w_state] <- 'W'
  list(classes = data.frame(date = df3$date, class = lab[vit],
                            stringsAsFactors = FALSE),
       loglik = -best$mod$minimum,
       means = paste0('step=', paste(round(step_m), collapse = '/'),
                      ';dist=', paste(round(dist_m), collapse = '/')))
}


# =============================================================================
# 4. WORKER (one track, all methods) ----
# =============================================================================

baselines_one <- function(track, info, p, root, cache_file, methods) {
  suppressPackageStartupMessages({
    for (pk in c('tidyverse', 'lubridate', 'sf', 'adehabitatLT', 'segclust2d')) {
      library(pk, character.only = TRUE)
    }
  })
  env <- new.env()
  sys.source(file.path(root, 'functions', 'segmentation_classification_functions.R'),
             envir = env)
  sys.source(file.path(root, 'functions', 'phenology_dates_functions.R'),
             envir = env)
  err_msg <- function(e) paste0(conditionMessage(e), ' [in: ',
                                deparse(conditionCall(e))[1], ']')

  df3 <- tryCatch(suppressWarnings(suppressMessages(prepare_track(track, info, p, env))),
                  error = err_msg)
  if (is.character(df3)) {
    return(lapply(setNames(methods, methods), function(m) list(ok = FALSE, error = df3)))
  }

  out <- list()
  for (m in methods) {
    t0 <- Sys.time()
    res <- tryCatch(suppressWarnings(suppressMessages(switch(m,
      segclust_simple = list(classes = baseline_segclust_simple(df3, p, env, cache_file)),
      migrateR = baseline_migrater(df3, p),
      hmm = baseline_hmm(df3, p)))), error = err_msg)
    if (is.character(res)) { out[[m]] <- list(ok = FALSE, error = res); next }
    cls <- res$classes
    cls$class <- final_classes(cls$class, p$lmin, env)
    dts <- if (m == 'migrateR') res$dates else dates_from_classes(cls$date, cls$class)
    info_m <- data.frame(ID = info$geo_trip, method = m,
                         run_time_s = round(as.numeric(difftime(Sys.time(), t0,
                                                                units = 'secs')), 1),
                         model = if (is.null(res$model)) NA else res$model,
                         details = if (!is.null(res$aic)) res$aic else
                           if (!is.null(res$means)) res$means else NA)
    out[[m]] <- list(ok = TRUE,
                     classes = data.frame(ID = info$geo_trip, cls),
                     dates = data.frame(ID = info$geo_trip, as.data.frame(dts)),
                     info = info_m)
  }
  out
}
