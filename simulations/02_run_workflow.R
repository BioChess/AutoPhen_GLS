# =============================================================================
# SIMULATIONS 02: RUN THE WORKFLOW ON THE SIMULATED TRACKS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Runs scripts 01 (segmentation and classification) and 03
#              (phenological dates) on the simulated tracks, fully automated
#              (no supervised correction), with the default parameters of
#              parameters.R.
#                - Tracks with GLS-like error: data_type = 'GLS' (equinox
#                  handling on).
#                - Tracks without error (GPS-like): data_type = 'GPS'.
#              The same runner as the sensitivity tests is used
#              (analyses/sensitivity/sensitivity_functions.R), so tracks can be
#              processed in parallel R sessions.
#
# Input:  simulations/data/sim_positions.csv, sim_metadata.csv
#         (from 01_simulate_tracks.R)
# Output: output/simulations/<gls|none>/ (Segmentation_results.csv,
#         Summary_segmentation.csv, Phenology_tracks.csv, Phenology.csv,
#         Failed_tracks.csv if any track failed, variant_info.csv)
#
# Run from the project root. About 40-60 s per track in one R session.
# =============================================================================

pacman::p_load('tidyverse', 'furrr')
source('analyses/sensitivity/sensitivity_functions.R')

# Parallel R sessions (1 = no parallel processing)
n_workers <- 4
# Run again groups that are already done?
overwrite <- FALSE
# Run again only the tracks that failed in a previous run (listed in
# output/simulations/<group>/Failed_tracks.csv) and add them to the results
rerun_failed <- FALSE

metadata <- read_csv(file.path('simulations', 'data', 'sim_metadata.csv'),
                     show_col_types = FALSE, col_types = cols(geo_trip = 'c'))

# Two groups of tracks with their settings (same structure as the sensitivity
# variants: changes to parameters.R)
groups <- list(
  gls  = list(group = 'Simulations', parameter = 'Observation error',
              value = 'GLS-like', default = '-', step = '01',
              changes = list(data_type = 'GLS',
                             positions_file = file.path('simulations', 'data',
                                                        'sim_positions.csv'))),
  none = list(group = 'Simulations', parameter = 'Observation error',
              value = 'None (GPS-like)', default = '-', step = '01',
              changes = list(data_type = 'GPS',
                             positions_file = file.path('simulations', 'data',
                                                        'sim_positions.csv')))
)
ids <- list(gls  = metadata$geo_trip[metadata$error %in% c('gls', 'gls_eq')],
            none = metadata$geo_trip[metadata$error == 'none'])

if (n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)
options(future.globals.maxSize = 2 * 1024^3)

out_root <- file.path('output', 'simulations')
for (g in names(groups)) {
  if (!rerun_failed) {
    run_variant(g, groups, ids = ids[[g]], metadata = metadata,
                n_workers = n_workers, overwrite = overwrite,
                out_root = out_root,
                cache_root = file.path(out_root, '_segclust_cache'))
    next
  }

  # --- Only the failed tracks ----------------------------------------------------
  failed_file <- file.path(out_root, g, 'Failed_tracks.csv')
  if (!file.exists(failed_file)) next
  failed <- read_csv(failed_file, show_col_types = FALSE, col_types = cols(ID = 'c'))
  cat('\nGroup', g, ':', nrow(failed), 'failed tracks to run again\n')
  tmp_root <- file.path(out_root, '_rerun')
  unlink(file.path(tmp_root, g), recursive = TRUE)
  run_variant(g, groups, ids = unique(failed$ID), metadata = metadata,
              n_workers = n_workers, overwrite = TRUE, out_root = tmp_root,
              cache_root = file.path(out_root, '_segclust_cache'))
  # Add the new results to the files of the group
  for (f in c('Segmentation_results.csv', 'Summary_segmentation.csv',
              'Phenology_tracks.csv', 'Phenology.csv')) {
    new_f <- file.path(tmp_root, g, f)
    if (!file.exists(new_f)) next
    old <- read_csv(file.path(out_root, g, f), show_col_types = FALSE,
                    col_types = cols(ID = 'c', .default = 'c'))
    add <- read_csv(new_f, show_col_types = FALSE, col_types = cols(ID = 'c', .default = 'c'))
    write_csv(bind_rows(old %>% filter(!ID %in% add$ID), add), file.path(out_root, g, f))
  }
  still <- file.path(tmp_root, g, 'Failed_tracks.csv')
  if (file.exists(still)) file.copy(still, failed_file, overwrite = TRUE) else unlink(failed_file)
  cat('Group', g, ':', ifelse(file.exists(still), nrow(read_csv(still, show_col_types = FALSE)), 0),
      'tracks still failing\n')
}
plan(sequential)
