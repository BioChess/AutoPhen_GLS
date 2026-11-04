# =============================================================================
# SENSITIVITY TESTS: FUNCTIONS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Functions to run the workflow (steps 01 and 03) with the
#              parameters of one variant, and to compare the results with the
#              reference (supervised) classification and dates.
#              Used by 02_run_variants.R and 03_evaluate_variants.R.
# =============================================================================


# =============================================================================
# 1. PARAMETERS OF A VARIANT ----
# =============================================================================

# Values derived from other parameters (as at the end of parameters.R)
derive_parameters <- function(p) {
  p$use_equinox      <- p$data_type == 'GLS'
  p$lmin             <- max(5, ceiling(p$min_segment_days * p$positions_per_day))
  p$min_stage_length <- max(1, ceiling(p$min_stage_days * p$positions_per_day))
  p
}

# Parameters of parameters.R, with the changes of one variant
variant_parameters <- function(v) {
  e <- new.env()
  sys.source('parameters.R', envir = e)
  p <- as.list(e)
  p$positions_file <- file.path('data', 'All_positions.csv')
  derive_parameters(modifyList(p, v$changes))
}

# Name of the segclust2d result saved for one track. The name contains all
# the settings that change the segmentation, so a saved result is only reused
# by variants with the same positions and the same segmentation settings
# (e.g. the variants of the classification thresholds).
segclust_key <- function(p) {
  paste0(tools::file_path_sans_ext(basename(p$positions_file)),
         '_step', p$adv$max_step,
         '_eq', ifelse(p$use_equinox, p$eq_window_days, 'no'),
         '_rep', p$adv$max_repeated, '-', p$adv$edge_positions,
         '_K', p$Kmax, '_lmin', p$lmin,
         '_ncl', min(p$ncluster), '-', max(p$ncluster))
}


# =============================================================================
# 2. WORKERS (one track; they are run in parallel R sessions) ----
# =============================================================================
# Each worker loads the packages and the functions of the workflow itself, so
# the same code works with one or several R sessions.

classify_one <- function(track, info, p, root, cache_file) {
  suppressPackageStartupMessages({
    for (pk in c('tidyverse', 'lubridate', 'sf', 'adehabitatLT', 'adehabitatHR',
                 'segclust2d', 'trajr', 'recipes', 'units')) {
      library(pk, character.only = TRUE)
    }
  })
  env <- new.env()
  sys.source(file.path(root, 'functions', 'segmentation_classification_functions.R'),
             envir = env)
  t0 <- Sys.time()
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      env$classify_track(track = track, info = info,
                         Kmax = p$Kmax, lmin = p$lmin, ncluster = p$ncluster,
                         eq_window_days = p$eq_window_days, adv = p$adv,
                         use_equinox = p$use_equinox, make_plots = FALSE,
                         segclust_cache = cache_file))),
    error = function(e) conditionMessage(e))
  if (is.character(res)) return(list(ok = FALSE, error = res))
  res$summary$run_time_s <- round(as.numeric(difftime(Sys.time(), t0,
                                                      units = 'secs')))
  res$summary$segclust_from_cache <- !is.null(cache_file) &&
    file.exists(cache_file) && file.mtime(cache_file) < t0
  res$plot <- NULL
  c(list(ok = TRUE), res)
}

phenology_one <- function(track.ph, info, p, root) {
  suppressPackageStartupMessages({
    for (pk in c('tidyverse', 'lubridate', 'sf', 'adehabitatHR', 'track2KBA')) {
      library(pk, character.only = TRUE)
    }
  })
  sf::sf_use_s2(FALSE)
  env <- new.env()
  sys.source(file.path(root, 'functions', 'phenology_dates_functions.R'),
             envir = env)
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      env$phenology_track(track.ph = track.ph, info = info,
                          min_stage_length = p$min_stage_length,
                          kde_level = p$kde_level, adv = p$adv,
                          use_equinox = p$use_equinox, world = NULL,
                          make_plot = FALSE))),
    error = function(e) conditionMessage(e))
  if (is.character(res)) return(list(ok = FALSE, error = res))
  list(ok = TRUE, dates = res$dates)
}

# Runs 'fun' for each element of 'args' (a list of argument lists), in
# parallel if n_workers > 1. Shows the progress every 'chunk' tracks.
run_tracks <- function(fun, args, n_workers, label) {
  chunk  <- max(1, 2 * n_workers)
  groups <- split(seq_along(args), ceiling(seq_along(args) / chunk))
  out <- vector('list', length(args))
  t0 <- Sys.time()
  for (g in groups) {
    out[g] <- furrr::future_map(args[g], function(a) do.call(fun, a),
                                .options = furrr::furrr_options(seed = 1L))
    cat(sprintf('  %s: %d/%d tracks (%.1f min)\n', label, max(g), length(args),
                as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
  }
  names(out) <- names(args)
  out
}


# =============================================================================
# 3. RUN ONE VARIANT ----
# =============================================================================

run_variant <- function(name, variants, ids, metadata, n_workers = 1,
                        out_root = file.path('output', 'sensitivity'),
                        cache_root = file.path('output', 'sensitivity', '_segclust_cache'),
                        overwrite = FALSE) {

  v <- variants[[name]]
  if (is.null(v)) stop('Variant "', name, '" not found in sensitivity_variants.R')
  p <- variant_parameters(v)
  out_dir <- file.path(out_root, name)
  root <- normalizePath('.')

  if (file.exists(file.path(out_dir, 'Phenology.csv')) && !overwrite) {
    cat('Variant', name, 'already done (overwrite = FALSE)\n')
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat('\n=== Variant', name, '|', v$parameter, '=', v$value,
      '(default', v$default, ') ===\n')
  t_start <- Sys.time()

  failed <- data.frame(ID = character(), step = character(), error = character())

  # --- Step 01 ------------------------------------------------------------------
  if (v$step == '01') {
    if (!file.exists(p$positions_file)) {
      stop('File not found: ', p$positions_file,
           '\nRun analyses/sensitivity/01_prepare_cdays_positions.R first')
    }
    positions <- readr::read_csv(p$positions_file, show_col_types = FALSE,
                                 col_types = readr::cols(geo_trip = 'c')) %>%
      dplyr::filter(geo_trip %in% ids)
    cache_dir <- file.path(cache_root, segclust_key(p))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    args <- lapply(ids, function(id) {
      list(track = positions[positions$geo_trip == id, ],
           info = metadata[metadata$geo_trip == id, ],
           p = p, root = root,
           cache_file = normalizePath(file.path(cache_dir, paste0(id, '.rds')),
                                      mustWork = FALSE))
    })
    names(args) <- ids
    res <- run_tracks(classify_one, args, n_workers, 'step 01')

    ok <- sapply(res, `[[`, 'ok')
    if (any(!ok)) {
      failed <- rbind(failed, data.frame(ID = ids[!ok], step = '01',
                                         error = sapply(res[!ok], `[[`, 'error')))
    }
    fix_id <- function(d) dplyr::mutate(d, ID = as.character(ID))
    readr::write_csv(dplyr::bind_rows(lapply(res[ok], function(r) fix_id(r$segments))),
                     file.path(out_dir, 'Segmentation_results.csv'))
    readr::write_csv(dplyr::bind_rows(lapply(res[ok], function(r) fix_id(r$summary))),
                     file.path(out_dir, 'Summary_segmentation.csv'))
    readr::write_csv(dplyr::bind_rows(lapply(res[ok], function(r) fix_id(r$track))),
                     file.path(out_dir, 'Phenology_tracks.csv'))
    tracks_file <- file.path(out_dir, 'Phenology_tracks.csv')
  } else {
    # Only step 03: classification of the main analysis
    tracks_file <- file.path(out_root, 'default', 'Phenology_tracks.csv')
    if (!file.exists(tracks_file)) stop('Run the variant "default" first')
  }

  # --- Step 03 ------------------------------------------------------------------
  tracks.ph <- readr::read_csv(tracks_file, show_col_types = FALSE,
                               col_types = readr::cols(ID = 'c'))
  ids3 <- intersect(ids, unique(tracks.ph$ID))
  args <- lapply(ids3, function(id) {
    list(track.ph = tracks.ph[tracks.ph$ID == id, ],
         info = metadata[metadata$geo_trip == id, ], p = p, root = root)
  })
  names(args) <- ids3
  res <- run_tracks(phenology_one, args, n_workers, 'step 03')
  ok <- sapply(res, `[[`, 'ok')
  if (any(!ok)) {
    failed <- rbind(failed, data.frame(ID = ids3[!ok], step = '03',
                                       error = sapply(res[!ok], `[[`, 'error')))
  }
  readr::write_csv(dplyr::bind_rows(lapply(res[ok], `[[`, 'dates')),
                   file.path(out_dir, 'Phenology.csv'))

  # --- Information of the run ---------------------------------------------------
  if (nrow(failed) > 0) readr::write_csv(failed, file.path(out_dir, 'Failed_tracks.csv'))
  info_run <- data.frame(variant = name, group = v$group, parameter = v$parameter,
                         value = v$value, default = v$default, step = v$step,
                         n_tracks = length(ids), n_failed = nrow(failed),
                         n_workers = n_workers,
                         total_time_min = round(as.numeric(difftime(
                           Sys.time(), t_start, units = 'mins')), 1),
                         date = format(Sys.time(), '%Y-%m-%d %H:%M'),
                         R = R.version.string)
  readr::write_csv(info_run, file.path(out_dir, 'variant_info.csv'))
  saveRDS(p, file.path(out_dir, 'parameters.rds'))
  cat('Variant', name, 'finished in', info_run$total_time_min, 'min;',
      nrow(failed), 'failed\n')
  invisible(info_run)
}


# =============================================================================
# 4. COMPARISON WITH THE REFERENCE ----
# =============================================================================

# Breeding / migration / wintering (mig1 and mig2 are both 'mig')
three_classes <- function(x) ifelse(x %in% c('B', 'W'), x, 'mig')

# Last departure from a wintering area (the tracks can have several W stages)
last_W_dep <- function(d) {
  cols <- grep('^W[0-9]+_dep', names(d), value = TRUE)
  if (length(cols) == 0) return(as.Date(rep(NA, nrow(d))))
  x <- as.data.frame(lapply(d[cols], function(z) as.Date(z)))
  as.Date(apply(x, 1, function(z) { z <- z[!is.na(z)]
    if (length(z)) as.character(max(as.Date(z))) else NA }))
}

# Four main dates of each track, with standard names
main_dates <- function(d, suffix = '') {
  g <- function(n) {
    col <- paste0(n, suffix)
    if (col %in% names(d)) as.Date(d[[col]]) else as.Date(rep(NA, nrow(d)))
  }
  w <- d; names(w) <- sub(paste0(suffix, '$'), '', names(w))
  data.frame(ID = as.character(d$ID), B1_dep = g('B1_dep'), W_arr = g('W1_arr'),
             W_dep = last_W_dep(w), B2_arr = g('B2_arr'))
}

date_names <- c('B1_dep', 'W_arr', 'W_dep', 'B2_arr')

# Per-track results of one variant
evaluate_variant <- function(out_dir, ref_class, ref_dates, class_dir = out_dir) {

  # Class of each position (for step '03' variants: the class of 'default')
  pt <- readr::read_csv(file.path(class_dir, 'Phenology_tracks.csv'),
                        show_col_types = FALSE, col_types = readr::cols(ID = 'c')) %>%
    dplyr::transmute(ID, date, new = three_classes(Breed2))
  cls <- dplyr::inner_join(ref_class, pt, by = c('ID', 'date')) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(n_positions = dplyr::n(), class_agreement = mean(ref == new),
                     .groups = 'drop')

  # Segments flagged for review, number of segments and run time
  seg <- readr::read_csv(file.path(class_dir, 'Segmentation_results.csv'),
                         show_col_types = FALSE, col_types = readr::cols(ID = 'c')) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(n_segments = dplyr::n(),
                     prop_check = mean(Check == 'Check'), .groups = 'drop')
  sm <- readr::read_csv(file.path(class_dir, 'Summary_segmentation.csv'),
                        show_col_types = FALSE, col_types = readr::cols(ID = 'c')) %>%
    dplyr::select(ID, run_time_s)

  # Dates
  ph <- readr::read_csv(file.path(out_dir, 'Phenology.csv'), show_col_types = FALSE,
                        col_types = readr::cols(ID = 'c'))
  new <- main_dates(ph, suffix = '_KDE')
  new$dates_check <- ph$Check
  dts <- dplyr::inner_join(ref_dates, new, by = 'ID', suffix = c('_ref', '_new'))
  for (n in date_names) {
    dts[[paste0('err_', n)]] <- as.numeric(dts[[paste0(n, '_new')]] -
                                             dts[[paste0(n, '_ref')]])
  }
  err <- as.matrix(dts[paste0('err_', date_names)])
  dts$mean_abs_error <- rowMeans(abs(err), na.rm = TRUE)
  dts$prop_within3   <- rowMeans(abs(err) <= 3, na.rm = TRUE)
  # Dates in the reference but missing in the variant
  has_ref <- !is.na(as.matrix(dts[paste0(date_names, '_ref')]))
  dts$n_missing <- rowSums(has_ref & is.na(err))

  cls %>%
    dplyr::full_join(seg, by = 'ID') %>%
    dplyr::full_join(sm, by = 'ID') %>%
    dplyr::full_join(dts %>% dplyr::select(ID, dates_check,
                                           dplyr::starts_with('err_'),
                                           mean_abs_error, prop_within3, n_missing),
                     by = 'ID') %>%
    dplyr::filter(ID %in% ref_class$ID)
}

# Bootstrap 95% CI of the mean of the paired differences (variant - default)
boot_ci <- function(x, n_boot = 2000) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(c(NA, NA))
  m <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  unname(quantile(m, c(0.025, 0.975)))
}
