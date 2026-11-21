# =============================================================================
# SIMULATIONS 03: COMPARISON WITH THE TRUE CLASSES AND DATES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Compares the automated output of the workflow with the truth of
#              the simulations.
#   Per track:   class agreement (proportion of positions with the true class:
#                breeding, migration or wintering; positions of a stopover are
#                correct as migration or as wintering),
#                signed error of each date (days; estimated - true; negative =
#                too early), failed tracks.
#   Per scenario (distance x movement x error) and date: number of dates,
#                bias (median signed error), median and mean absolute error,
#                90th percentile and maximum absolute error, percentage within
#                1, 3 and 7 days, missing dates.
#   Per segment: agreement score vs proportion of positions correctly
#                classified (does the agreement score predict errors?).
#
# Input:  simulations/data/ (truth), output/simulations/<gls|none>/
# Output: output/simulations/_results/ (tables; figures in 04_simulation_figures.R)
# =============================================================================

pacman::p_load('tidyverse')
source('analyses/sensitivity/sensitivity_functions.R')   # main_dates(), three_classes()

sim_dir <- file.path('simulations', 'data')
out_dir <- file.path('output', 'simulations')
res_dir <- file.path(out_dir, '_results')
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

metadata    <- read_csv(file.path(sim_dir, 'sim_metadata.csv'), show_col_types = FALSE,
                        col_types = cols(geo_trip = 'c'))
truth_dates <- read_csv(file.path(sim_dir, 'sim_truth_dates.csv'), show_col_types = FALSE) %>%
  rename(ID = sim_id)
truth_pos   <- read_csv(file.path(sim_dir, 'sim_truth_positions.csv'), show_col_types = FALSE,
                        col_types = cols(geo_trip = 'c')) %>%
  transmute(ID = geo_trip, date = Date_Time, true = stage, leg)

groups <- c('gls', 'none')
groups <- groups[file.exists(file.path(out_dir, groups, 'Phenology.csv'))]
read_group <- function(file) {
  map_dfr(groups, ~ read_csv(file.path(out_dir, .x, file), show_col_types = FALSE,
                             col_types = cols(ID = 'c')))
}
scen <- metadata %>% select(ID = geo_trip, distance, movement, error) %>%
  left_join(truth_dates %>% select(ID, equinox_overlap), by = 'ID')
lev <- function(d) d %>% mutate(
  distance = factor(distance, levels = c('short', 'medium', 'long')),
  movement = factor(movement, levels = c('base', 'foraging', 'stopover', 'two_W')),
  error    = factor(error, levels = c('none', 'gls', 'gls_eq')))


# =============================================================================
# 1. CLASSIFICATION ----
# =============================================================================

pt <- read_group('Phenology_tracks.csv') %>% transmute(ID, date, est = three_classes(Breed2))
# A stopover classified as wintering is accepted: it is a short stationary
# period, and what matters is that the real wintering area is also detected
# (this is checked with the dates below)
cls <- inner_join(truth_pos, pt, by = c('ID', 'date')) %>%
  mutate(correct = true == est | (leg == 'stopover' & est == 'W')) %>%
  group_by(ID) %>%
  summarise(class_agreement = mean(correct), n_positions = n(),
            # stopover detected as a wintering stage (most of its positions W)
            stopover_as_W = if (any(leg == 'stopover')) mean(est[leg == 'stopover'] == 'W') > 0.5 else NA,
            .groups = 'drop')


# =============================================================================
# 2. DATES ----
# =============================================================================

ph  <- read_group('Phenology.csv')
est <- main_dates(ph, suffix = '_KDE')
tru <- truth_dates %>% transmute(ID, B1_dep, W_arr, W_dep, B2_arr)

# Arrival at the wintering area: the first estimated wintering stage that
# overlaps the true wintering period. A stopover detected as an earlier
# wintering stage is then not taken as the arrival (it is only an error if
# the real wintering area is not detected).
w_stages <- ph %>%
  select(ID, matches('^W[0-9]+_(arr|dep)_KDE$')) %>%
  pivot_longer(-ID, names_to = c('stage', '.value'),
               names_pattern = '^(W[0-9]+)_(arr|dep)_KDE$') %>%
  mutate(arr = as.Date(arr), dep = as.Date(dep)) %>%
  filter(!is.na(arr))
w_arr_overlap <- w_stages %>%
  inner_join(tru %>% select(ID, W_arr_t = W_arr, W_dep_t = W_dep), by = 'ID') %>%
  filter(arr <= W_dep_t, is.na(dep) | dep >= W_arr_t) %>%
  group_by(ID) %>% summarise(W_arr_ov = min(arr), .groups = 'drop')
est <- est %>% left_join(w_arr_overlap, by = 'ID') %>%
  mutate(W_arr = coalesce(W_arr_ov, W_arr)) %>% select(-W_arr_ov)
dts <- inner_join(tru, est, by = 'ID', suffix = c('_true', '_est'))
err <- map_dfr(date_names, function(n) {
  tibble(ID = dts$ID, date = n,
         diff = as.numeric(dts[[paste0(n, '_est')]] - dts[[paste0(n, '_true')]]))
})

# Failed tracks (any step)
failed <- map_dfr(groups, function(g) {
  f <- file.path(out_dir, g, 'Failed_tracks.csv')
  if (file.exists(f)) read_csv(f, show_col_types = FALSE, col_types = cols(ID = 'c')) else NULL
})
if (nrow(failed) == 0) failed <- tibble(ID = character(), step = character(), error = character())
write_csv(failed, file.path(res_dir, 'sim_failed_tracks.csv'))
cat(nrow(failed), 'failed tracks of', nrow(metadata), '\n')


# =============================================================================
# 3. TABLES ----
# =============================================================================

per_track <- scen %>%
  left_join(cls, by = 'ID') %>%
  left_join(err %>% pivot_wider(names_from = date, values_from = diff,
                                names_prefix = 'err_'), by = 'ID') %>%
  mutate(failed = ID %in% failed$ID) %>% lev()
write_csv(per_track, file.path(res_dir, 'sim_per_track.csv'))

summ_dates <- function(d) d %>% summarise(
  n_dates      = sum(!is.na(diff)),
  missing      = sum(is.na(diff)),
  bias_median  = median(diff, na.rm = TRUE),
  abs_median   = median(abs(diff), na.rm = TRUE),
  abs_mean     = round(mean(abs(diff), na.rm = TRUE), 1),
  abs_p90      = unname(quantile(abs(diff), 0.9, na.rm = TRUE)),
  abs_max      = max(abs(diff), na.rm = TRUE),
  within1_pct  = round(100 * mean(abs(diff) <= 1, na.rm = TRUE), 1),
  within3_pct  = round(100 * mean(abs(diff) <= 3, na.rm = TRUE), 1),
  within7_pct  = round(100 * mean(abs(diff) <= 7, na.rm = TRUE), 1),
  .groups = 'drop')

err_s <- err %>% left_join(scen, by = 'ID') %>% lev() %>%
  mutate(date = factor(date, levels = date_names))

# By scenario and date
tab_dates <- err_s %>% group_by(error, distance, movement, date) %>% summ_dates()
write_csv(tab_dates, file.path(res_dir, 'sim_dates_by_scenario.csv'))
# By error and date (summary for the main text)
tab_dates_error <- err_s %>% group_by(error, date) %>% summ_dates()
# GLS tracks by the equinox overlapped by a migration
write_csv(err_s %>% filter(error != 'none') %>%
            group_by(error, equinox_overlap, date) %>% summ_dates(),
          file.path(res_dir, 'sim_dates_by_equinox.csv'))
# Stopovers detected as wintering stages
write_csv(per_track %>% filter(movement == 'stopover') %>%
            group_by(error, distance) %>%
            summarise(n_tracks = n(), stopover_as_W_pct = round(100 * mean(stopover_as_W, na.rm = TRUE), 1),
                      W_arr_within3_pct = round(100 * mean(abs(err_W_arr) <= 3, na.rm = TRUE), 1),
                      .groups = 'drop'),
          file.path(res_dir, 'sim_stopovers.csv'))
write_csv(tab_dates_error, file.path(res_dir, 'sim_dates_by_error.csv'))
print(tab_dates_error, n = Inf, width = Inf)

# Classification by scenario
tab_class <- per_track %>% group_by(error, distance, movement) %>%
  summarise(n_tracks = n(), failed = sum(failed),
            class_agreement_pct = round(100 * weighted.mean(class_agreement, n_positions,
                                                            na.rm = TRUE), 1),
            class_median_track = round(100 * median(class_agreement, na.rm = TRUE), 1),
            .groups = 'drop')
write_csv(tab_class, file.path(res_dir, 'sim_class_by_scenario.csv'))
print(per_track %>% group_by(error) %>%
        summarise(class_agreement_pct = round(100 * weighted.mean(class_agreement, n_positions,
                                                                  na.rm = TRUE), 1),
                  failed = sum(failed)))


# =============================================================================
# 4. AGREEMENT SCORE AND SEGMENTS FLAGGED FOR REVIEW ----
# =============================================================================
# For each segment: proportion of its positions with the true class. A segment
# is 'wrong' if most of its positions are misclassified.

seg <- read_group('Segmentation_results.csv')
# Positions of each segment: rows begin to end of the track (same rows as
# Phenology_tracks.csv)
pos <- read_group('Phenology_tracks.csv') %>%
  group_by(ID) %>% mutate(row = row_number()) %>% ungroup() %>%
  select(ID, date, row) %>%
  inner_join(truth_pos, by = c('ID', 'date'))
seg_true <- seg %>%
  select(ID, begin, end, Breed, Agreement, Check, Check_reason) %>%
  mutate(seg_n = row_number()) %>%
  inner_join(pos, by = join_by(ID, begin <= row, end >= row)) %>%
  group_by(seg_n, ID, Breed, Agreement, Check, Check_reason) %>%
  summarise(prop_correct = mean(three_classes(Breed) == true |
                                (leg == 'stopover' & three_classes(Breed) == 'W')),
            .groups = 'drop') %>%
  mutate(wrong = prop_correct < 0.5) %>%
  left_join(scen, by = 'ID') %>% lev()
write_csv(seg_true, file.path(res_dir, 'sim_segments.csv'))

# Agreement score vs wrong segments, and flagged segments
tab_agree <- seg_true %>%
  mutate(agreement_class = cut(Agreement, c(-0.01, 0.25, 0.5, 0.75, 1),
                               labels = c('0-0.25', '0.25-0.5', '0.5-0.75', '0.75-1'))) %>%
  group_by(error, agreement_class) %>%
  summarise(n_segments = n(), wrong_pct = round(100 * mean(wrong), 1), .groups = 'drop')
write_csv(tab_agree, file.path(res_dir, 'sim_agreement_vs_error.csv'))
tab_check <- seg_true %>% group_by(error) %>%
  summarise(n_segments = n(),
            flagged_pct = round(100 * mean(Check == 'Check'), 1),
            wrong_pct = round(100 * mean(wrong), 1),
            wrong_flagged_pct = round(100 * mean(Check[wrong] == 'Check'), 1),
            wrong_if_flagged_pct = round(100 * mean(wrong[Check == 'Check']), 1),
            wrong_if_ok_pct = round(100 * mean(wrong[Check != 'Check']), 1),
            .groups = 'drop')
write_csv(tab_check, file.path(res_dir, 'sim_flagged_segments.csv'))
print(tab_check, width = Inf)


# Figures: simulations/04_simulation_figures.R
