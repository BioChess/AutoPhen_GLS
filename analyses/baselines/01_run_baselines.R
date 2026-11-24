# =============================================================================
# BASELINE METHODS 01: RUN THE BASELINES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Runs three simple alternatives to the workflow (see
#              baseline_functions.R) on the simulated tracks and on the real
#              tracks with a supervised reference (Calonectris), with the same
#              positions and cleaning as the workflow:
#                segclust_simple  segclust2d + simple rule (NSCD limit and
#                                 net speed)
#                migrateR         NSD models (migrant / mixmig dates)
#                hmm              3-state hidden Markov model (momentuHMM)
#              The segclust2d results of the workflow are read from its cache
#              (output/simulations/_segclust_cache and
#              output/sensitivity/_segclust_cache), so the segmentation is not
#              calculated again.
#
# Input:  simulations/data/ (sim_positions.csv, sim_metadata.csv)
#         data/All_positions.csv, data/All_metadata.csv,
#         data/reference/CAL_tracks_supervised.csv (IDs of the real tracks)
# Output: output/baselines/<dataset>/<method>_classes.csv  (ID, date, class)
#         output/baselines/<dataset>/<method>_dates.csv    (ID, B1_dep, W_arr,
#                                                           W_dep, B2_arr)
#         output/baselines/<dataset>/<method>_info.csv     (run time, model)
#         output/baselines/<dataset>/Failed_tracks.csv     (if any)
#         <dataset> = sim_gls (GLS-like error), sim_none (no error), real
#
# Run from the project root. Needs momentuHMM (CRAN) and migrateR
# (remotes::install_github('dbspitz/migrateR/migrateR')).
# =============================================================================

pacman::p_load('tidyverse', 'furrr')
source('analyses/sensitivity/sensitivity_functions.R')   # parameters, runner
source('analyses/baselines/baseline_functions.R')


# =============================================================================
# 1. SETTINGS ----
# =============================================================================

# Datasets and methods to run
datasets <- c('sim_gls', 'sim_none', 'real')
methods  <- c('segclust_simple', 'migrateR', 'hmm')

# Parallel R sessions (1 = no parallel processing)
n_workers <- 4

# Run again datasets that are already done?
overwrite <- FALSE

# Settings of the baselines (the rest are the defaults of parameters.R)
baseline_settings <- list(
  # segclust_simple: a non-breeding segment is migration if its net speed is
  # higher than this (km/day)
  simple_speed_km_day = 100,
  # migrateR: starting values of the time of the outward movement (theta) and
  # of the stay (rho), as proportions of the track duration (all combinations
  # are fitted; the best fit of each model is kept); maximum duration
  # parameter of each movement (phi, days); dates at 5% and 95% of the change
  migrater_theta_start = c(0.2, 0.3, 0.4, 0.5),
  migrater_rho_start   = c(0.3, 0.45, 0.6),
  migrater_phi_max     = 21,
  migrater_p           = 0.05,
  # hmm: number of starting values (the first one from quantiles of the data,
  # the others at random around it)
  hmm_n_starts = 3)


# =============================================================================
# 2. TRACKS OF EACH DATASET ----
# =============================================================================

sim_meta <- function() {
  read_csv(file.path('simulations', 'data', 'sim_metadata.csv'),
           show_col_types = FALSE, col_types = cols(geo_trip = 'c'))
}
dataset_setup <- function(d) {
  if (d %in% c('sim_gls', 'sim_none')) {
    md <- sim_meta()
    err <- if (d == 'sim_gls') c('gls', 'gls_eq') else 'none'
    list(metadata = md,
         ids = md$geo_trip[md$error %in% err],
         changes = list(data_type = if (d == 'sim_gls') 'GLS' else 'GPS',
                        positions_file = file.path('simulations', 'data',
                                                   'sim_positions.csv')),
         cache_root = file.path('output', 'simulations', '_segclust_cache'))
  } else {
    md <- read_csv(file.path('data', 'All_metadata.csv'), show_col_types = FALSE,
                   col_types = cols(geo_trip = 'c'))
    ids <- read_csv(file.path('data', 'reference', 'CAL_tracks_supervised.csv'),
                    show_col_types = FALSE, col_select = 'ID',
                    col_types = cols(ID = 'c')) %>% distinct(ID) %>% pull(ID)
    list(metadata = md, ids = intersect(ids, md$geo_trip),
         changes = list(data_type = 'GLS',
                        positions_file = file.path('data', 'All_positions.csv')),
         cache_root = file.path('output', 'sensitivity', '_segclust_cache'))
  }
}


# =============================================================================
# 3. RUN ----
# =============================================================================

if (n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)
options(future.globals.maxSize = 2 * 1024^3)
root <- normalizePath('.')

for (d in datasets) {
  out_dir <- file.path('output', 'baselines', d)
  if (all(file.exists(file.path(out_dir, paste0(methods, '_dates.csv')))) && !overwrite) {
    cat('Dataset', d, 'already done (overwrite = FALSE)\n'); next
  }
  s <- dataset_setup(d)
  if (!file.exists(s$changes$positions_file)) {
    cat('Dataset', d, ': file not found,', s$changes$positions_file, '\n'); next
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p <- modifyList(variant_parameters(list(changes = s$changes)), baseline_settings)
  cache_dir <- file.path(s$cache_root, segclust_key(p))
  if (!dir.exists(cache_dir)) {
    cat('  No segclust2d cache for', d, '(', cache_dir, '): segmentation will be',
        'calculated again (slower)\n')
  }
  cat('\n=== Dataset', d, ':', length(s$ids), 'tracks ===\n')
  t0 <- Sys.time()

  positions <- read_csv(s$changes$positions_file, show_col_types = FALSE,
                        col_types = cols(geo_trip = 'c')) %>%
    filter(geo_trip %in% s$ids)
  args <- lapply(s$ids, function(id) {
    list(track = positions[positions$geo_trip == id, ],
         info = s$metadata[s$metadata$geo_trip == id, ],
         p = p, root = root,
         cache_file = normalizePath(file.path(cache_dir, paste0(id, '.rds')),
                                    mustWork = FALSE),
         methods = methods)
  })
  names(args) <- s$ids
  res <- run_tracks(baselines_one, args, n_workers, d)

  failed <- list()
  for (m in methods) {
    r  <- lapply(res, `[[`, m)
    ok <- sapply(r, `[[`, 'ok')
    if (any(!ok)) failed[[m]] <- data.frame(ID = names(r)[!ok], method = m,
                                            error = sapply(r[!ok], `[[`, 'error'))
    write_csv(bind_rows(lapply(r[ok], `[[`, 'classes')),
              file.path(out_dir, paste0(m, '_classes.csv')))
    write_csv(bind_rows(lapply(r[ok], `[[`, 'dates')),
              file.path(out_dir, paste0(m, '_dates.csv')))
    write_csv(bind_rows(lapply(r[ok], `[[`, 'info')),
              file.path(out_dir, paste0(m, '_info.csv')))
    cat(' ', m, ':', sum(ok), 'tracks done,', sum(!ok), 'failed\n')
  }
  failed <- bind_rows(failed)
  if (nrow(failed) > 0) write_csv(failed, file.path(out_dir, 'Failed_tracks.csv'))
  saveRDS(p, file.path(out_dir, 'parameters.rds'))
  cat('Dataset', d, 'finished in',
      round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), 'min\n')
}
