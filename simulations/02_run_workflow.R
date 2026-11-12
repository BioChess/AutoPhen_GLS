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

for (g in names(groups)) {
  run_variant(g, groups, ids = ids[[g]], metadata = metadata,
              n_workers = n_workers, overwrite = overwrite,
              out_root = file.path('output', 'simulations'),
              cache_root = file.path('output', 'simulations', '_segclust_cache'))
}
plan(sequential)
