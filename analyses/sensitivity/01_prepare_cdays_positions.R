# =============================================================================
# SENSITIVITY TESTS 01: POSITIONS WITH OTHER TREATMENTS OF THE COLONY DAYS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Runs analyses/prepare_full_dataset.R once for each treatment of
#              the colony days of Calonectris (kept, removed, anchored at the
#              colony, random point within 50 or 200 km, and random points
#              within 100 km with other seeds). The main analysis
#              uses random points within 100 km (data/All_positions.csv,
#              which is not changed).
#
# Input:  the same files as analyses/prepare_full_dataset.R
# Output: data/sensitivity/All_positions_cdays_<name>.csv
#         data/sensitivity/checks_<name>/ (cleaning summary of each version)
#
# Run from the project root. Only needed for the 'cdays_' variants.
# =============================================================================

source('analyses/sensitivity/sensitivity_variants.R')
dir.create(file.path('data', 'sensitivity'), recursive = TRUE, showWarnings = FALSE)

# Versions to prepare (names(cdays_variants) = all). Versions already
# prepared are skipped.
to_prepare <- names(cdays_variants)

for (name in to_prepare) {
  if (file.exists(cdays_file(name))) {
    cat('\n', cdays_file(name), 'already exists\n'); next
  }
  cat('\n=== Colony days:', name, '===\n')
  sens_cdays <- c(cdays_variants[[name]],
                  list(positions_file = cdays_file(name),
                       check_dir = file.path('data', 'sensitivity',
                                             paste0('checks_', name))))
  # local environment: the objects of the script do not stay in the session
  source('analyses/prepare_full_dataset.R', local = new.env())
}
if (exists("sens_cdays")) rm(sens_cdays)
