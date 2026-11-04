# =============================================================================
# SENSITIVITY TESTS 02: RUN THE VARIANTS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Runs the workflow (steps 01 and 03) for the variants listed in
#              sensitivity_variants.R, only for the tracks of the reference
#              (Calonectris tracks with supervised classification).
#              Each variant is saved in output/sensitivity/<variant>/ with the
#              same files as scripts 01 and 03 (Phenology_tracks.csv,
#              Segmentation_results.csv, Summary_segmentation.csv,
#              Phenology.csv) plus variant_info.csv (run time, number of
#              failed tracks) and parameters.rds.
#
#              Run 'default' first: the variants that only change the dates
#              (step '03') use its classification. Variants that are already
#              done are skipped (unless overwrite <- TRUE).
#
# Input:  data/All_metadata.csv, data/All_positions.csv (and the files of the
#         colony-day variants, see 01_prepare_cdays_positions.R)
#         data/reference/CAL_tracks_supervised.csv (IDs of the tracks)
# Output: output/sensitivity/<variant>/
#
# Run from the project root (open AutoPhen_GLS.Rproj).
# =============================================================================

pacman::p_load('tidyverse', 'furrr')
source('analyses/sensitivity/sensitivity_variants.R')
source('analyses/sensitivity/sensitivity_functions.R')


# =============================================================================
# 1. SETTINGS ----
# =============================================================================

# Variants to run (names in sensitivity_variants.R). Examples:
#   'default'                       main analysis (run it first)
#   c('kmax_25', 'kmax_35')         several variants, one after the other
#   names(variants)                 all variants
# to_run <- c('default')
# to_run <- c('kde_75', 'kde_95', 'bw_0.5', 'bw_1.5', 'dep_first',
            # 'minstage_1.5', 'minstage_5')
to_run <- c('thr_sa_0.4', 'thr_sa_0.6', 'thr_sin_0.8', 'thr_sin_1.2',
            'r5_0.65', 'r5_0.85')
# Parallel R sessions (tracks are processed at the same time). With 6 cores
# and 16 GB of RAM, 4-5 is a good value. 1 = no parallel processing.
n_workers <- 4

# Run again variants that are already done?
overwrite <- FALSE

# Tracks of the reference
ref_class_file <- file.path('data', 'reference', 'CAL_tracks_supervised.csv')


# =============================================================================
# 2. RUN ----
# =============================================================================

metadata <- read_csv(file.path('data', 'All_metadata.csv'), show_col_types = FALSE,
                     col_types = cols(geo_trip = 'c'))
ids <- read_csv(ref_class_file, show_col_types = FALSE, col_select = 'ID',
                col_types = cols(ID = 'c')) %>%
  distinct(ID) %>% pull(ID)
ids <- intersect(ids, metadata$geo_trip)
cat(length(ids), 'tracks of the reference\n')

if (n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)
options(future.globals.maxSize = 2 * 1024^3)

for (name in to_run) {
  run_variant(name, variants, ids = ids, metadata = metadata,
              n_workers = n_workers, overwrite = overwrite)
}
plan(sequential)
