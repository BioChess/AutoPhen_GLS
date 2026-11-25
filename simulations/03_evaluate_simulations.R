# =============================================================================
# SIMULATIONS 03: COMPARISON WITH THE TRUE CLASSES AND DATES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Compares the automated output of the workflow with the truth of
#              the simulations.
#   Per track:   class agreement (proportion of positions with the true class:
#                breeding, migration or wintering; a stopover is a short
#                wintering stage, so its true class is wintering; the movement
#                between two wintering stages is migration, as in the
#                workflow),
#                signed error of each date (days; estimated - true; negative =
#                too early), failed tracks.
#   Per scenario (distance x movement x error) and date: number of dates,
#                bias (median signed error), median and mean absolute error,
#                90th percentile and maximum absolute error, percentage within
#                1, 3 and 7 days, missing dates. The wintering dates are those
#                of the total wintering period (arrival at the first wintering
#                area, departure from the last one).
#   Each wintering area (WA: stopover, W1, W2): % detected, errors of its
#                arrival and departure dates, and % of tracks with the true
#                number of WAs (sim_wintering_areas.csv). The main WA (the
#                longest true WA) errors are also in sim_per_track.csv.
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

# The segments are evaluated with the same true stage (true_seg = true).
truth_pos <- truth_pos %>% mutate(true_seg = true)

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
cls <- inner_join(truth_pos, pt, by = c('ID', 'date')) %>%
  group_by(ID) %>%
  summarise(class_agreement = mean(true == est), n_positions = n(),
            # stopover detected (most of its positions classified as W) and
            # its duration (days)
            stopover_detected = if (any(leg == 'stopover')) mean(est[leg == 'stopover'] == 'W') > 0.5 else NA,
            stopover_days = sum(leg == 'stopover') / 2,
            .groups = 'drop')


# =============================================================================
# 2. DATES ----
# =============================================================================

ph  <- read_group('Phenology.csv')
est <- main_dates(ph, suffix = '_KDE')
tru <- truth_dates %>% transmute(ID, B1_dep, W_arr, W_dep, B2_arr)

dts <- inner_join(tru, est, by = 'ID', suffix = c('_true', '_est'))
err <- map_dfr(date_names, function(n) {
  tibble(ID = dts$ID, date = n,
         diff = as.numeric(dts[[paste0(n, '_est')]] - dts[[paste0(n, '_true')]]))
})

# Each wintering area (WA). True WAs: runs of true wintering positions (the
# stopover, W1 and W2). Estimated WAs: wintering stages of the workflow
# (W1, W2, ... after the short-stage rule, as in step 03), with their KDE
# dates. Each true WA is matched with the estimated WA that shares most of
# its positions; 'detected' = a matched WA exists. The main WA is the longest
# true WA.
lmin <- readRDS(file.path(out_dir, groups[1], 'parameters.rds'))$lmin
check_phen_env <- new.env()
sys.source(file.path('functions', 'phenology_dates_functions.R'), envir = check_phen_env)
wa_runs <- function(d, col) {
  d %>% arrange(ID, date) %>% group_by(ID) %>%
    mutate(run = cumsum(.data[[col]] != lag(.data[[col]], default = first(.data[[col]])))) %>%
    filter(.data[[col]] == 'W') %>%
    mutate(k = dense_rank(run)) %>% ungroup() %>% select(ID, date, k)
}
true_wa <- truth_pos %>% wa_runs('true') %>% rename(k_true = k)
true_wa_dates <- truth_pos %>% arrange(ID, date) %>% group_by(ID) %>%
  mutate(next_date = lead(date)) %>% ungroup() %>%
  inner_join(true_wa, by = c('ID', 'date')) %>%
  group_by(ID, k_true) %>%
  summarise(arr_true = as.Date(min(date)),
            dep_true = as.Date(max(next_date, na.rm = TRUE)),
            days = as.numeric(difftime(max(date), min(date), units = 'days')),
            .groups = 'drop_last') %>%
  mutate(n_WA_true = n(), main = days == max(days)) %>% ungroup()
est_stage <- read_group('Phenology_tracks.csv') %>% select(ID, date, Breed2) %>%
  mutate(Breed2 = three_classes(Breed2)) %>% group_by(ID) %>% arrange(date, .by_group = TRUE) %>%
  group_modify(~ check_phen_env$check_phen(.x, 'Breed2', min_length = lmin)) %>% ungroup()
est_wa <- est_stage %>% wa_runs('Breed2') %>% rename(k_est = k)
match <- inner_join(true_wa, est_wa, by = c('ID', 'date')) %>%
  count(ID, k_true, k_est) %>% group_by(ID, k_true) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(ID, k_true, k_est)
est_wa_dates <- ph %>% select(ID, matches('^W[0-9]+_(arr|dep)_KDE$')) %>%
  pivot_longer(-ID, names_to = c('k_est', 'what'), names_pattern = 'W([0-9]+)_(arr|dep)_KDE',
               values_to = 'value') %>%
  mutate(k_est = as.integer(k_est), value = as.Date(value)) %>%
  pivot_wider(names_from = what, values_from = value)
n_wa_est <- est_wa %>% distinct(ID, k_est) %>% count(ID, name = 'n_WA_est')
wa_eval <- true_wa_dates %>% left_join(match, by = c('ID', 'k_true')) %>%
  left_join(est_wa_dates, by = c('ID', 'k_est')) %>%
  left_join(n_wa_est, by = 'ID') %>%
  mutate(n_WA_est = replace_na(n_WA_est, 0L), detected = !is.na(k_est),
         err_arr = as.numeric(arr - arr_true), err_dep = as.numeric(dep - dep_true))
main_wa <- wa_eval %>% filter(main) %>% group_by(ID) %>% slice(1) %>% ungroup() %>%
  transmute(ID, n_WA_true, n_WA_est, err_WA_main_arr = err_arr, err_WA_main_dep = err_dep)

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
  left_join(main_wa, by = 'ID') %>%
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

# Each wintering area: detection and dates, by scenario and type of WA
# (stopover; W1 and W2 of the tracks with two wintering areas; the only WA of
# the other tracks), and % of tracks with the true number of WAs
wa_tab <- wa_eval %>% left_join(scen, by = 'ID') %>%
  left_join(truth_pos %>% filter(true == 'W') %>% inner_join(true_wa, by = c('ID', 'date')) %>%
              distinct(ID, k_true, leg), by = c('ID', 'k_true')) %>%
  mutate(wa_type = case_when(leg == 'stopover' ~ 'stopover',
                             movement == 'two_W' & leg == 'W1' ~ 'first of two',
                             movement == 'two_W' & leg == 'W2' ~ 'second of two',
                             TRUE ~ 'main wintering area')) %>%
  lev()
write_csv(wa_tab, file.path(res_dir, 'sim_wintering_areas_per_track.csv'))
tab_wa <- wa_tab %>% pivot_longer(c(err_arr, err_dep), names_to = 'date', values_to = 'diff',
                                  names_prefix = 'err_') %>%
  group_by(error, movement, wa_type, date) %>%
  summ_dates() %>%
  left_join(wa_tab %>% group_by(error, movement, wa_type) %>%
              summarise(detected_pct = round(100 * mean(detected), 1), .groups = 'drop'),
            by = c('error', 'movement', 'wa_type')) %>%
  left_join(wa_tab %>% distinct(ID, error, movement, n_WA_true, n_WA_est) %>%
              group_by(error, movement) %>%
              summarise(n_WA_correct_pct = round(100 * mean(n_WA_est == n_WA_true), 1),
                        n_WA_more_pct = round(100 * mean(n_WA_est > n_WA_true), 1),
                        n_WA_fewer_pct = round(100 * mean(n_WA_est < n_WA_true), 1),
                        .groups = 'drop'), by = c('error', 'movement'))
write_csv(tab_wa, file.path(res_dir, 'sim_wintering_areas.csv'))

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
# Detection of short stationary stages (stopovers, 5-15 days) by duration:
# a limitation of the segmentation for short stages
write_csv(per_track %>% filter(movement == 'stopover') %>%
            mutate(stopover_duration = cut(stopover_days, c(0, 7.5, 10, 12.5, 16),
                                           labels = c('5-7.5', '7.5-10', '10-12.5', '12.5-15'))) %>%
            group_by(error, stopover_duration) %>%
            summarise(n_tracks = n(),
                      detected_pct = round(100 * mean(stopover_detected, na.rm = TRUE), 1),
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
  summarise(prop_correct = mean(three_classes(Breed) == true_seg), .groups = 'drop') %>%
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
