# =============================================================================
# SENSITIVITY TESTS 03: COMPARISON WITH THE REFERENCE
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Compares the results of each variant with the reference
#              (classification and dates after the supervised review) and
#              with the main analysis ('default').
#
#   For each track:
#     class_agreement  proportion of positions with the same class
#                      (breeding / migration / wintering) as the reference
#     err_<date>       difference (days) with the reference date
#                      (B1_dep, W_arr, W_dep = last wintering departure, B2_arr)
#     mean_abs_error   mean absolute difference of the four dates (days)
#     prop_within3     proportion of the dates within 3 days of the reference
#     prop_check       proportion of segments flagged for review
#     n_segments, run_time_s
#
#   For each variant (sens_summary.csv): the same values for all tracks, and
#   the mean paired difference with 'default' (variant - default, the same
#   tracks) with a bootstrap 95% confidence interval.
#
#   OPTIONAL: if the dates of a second analyst are available
#   (data/reference/Phen_other_analyst.xlsx), the dates of each variant are
#   also compared with them, and the difference between the two analysts is
#   calculated (variability between people).
#
# Input:  output/sensitivity/<variant>/
#         data/reference/CAL_tracks_supervised.csv (ID, date, Breed2)
#         data/reference/CAL_phen_supervised.csv   (ID, B1_dep, W1_arr, ...)
#         data/reference/Phen_other_analyst.xlsx   (optional)
# Output: output/sensitivity/_results/sens_per_track.csv
#         output/sensitivity/_results/sens_summary.csv
#         output/sensitivity/_results/sens_dates.csv
#         output/sensitivity/_results/other_analyst.csv (optional)
# =============================================================================

pacman::p_load('tidyverse', 'readxl')
source('analyses/sensitivity/sensitivity_variants.R')
source('analyses/sensitivity/sensitivity_functions.R')
set.seed(1)

out_root <- file.path('output', 'sensitivity')
res_dir  <- file.path(out_root, '_results')
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 1. REFERENCE ----
# =============================================================================

ref_class <- read_csv(file.path('data', 'reference', 'CAL_tracks_supervised.csv'),
                      show_col_types = FALSE,
                      col_types = cols(ID = 'c', date = col_datetime())) %>%
  transmute(ID, date, ref = three_classes(Breed2))

ref_dates <- read_csv(file.path('data', 'reference', 'CAL_phen_supervised.csv'),
                      show_col_types = FALSE, col_types = cols(ID = 'c')) %>%
  main_dates()

other_file <- file.path('data', 'reference', 'Phen_other_analyst.xlsx')
other <- NULL
if (file.exists(other_file)) {
  other_raw <- read_excel(other_file) %>%
    mutate(ID = as.character(geo_trip))
  other <- other_raw %>%
    group_by(Analyst) %>%
    group_modify(~ main_dates(.x)) %>%
    ungroup()
}


# =============================================================================
# 2. RESULTS OF EACH VARIANT ----
# =============================================================================

done <- names(variants)[file.exists(file.path(out_root, names(variants),
                                              'Phenology.csv'))]
if (!'default' %in% done) stop('The variant "default" is not done yet')
cat('Variants done:', paste(done, collapse = ', '), '\n')

per_track <- map_dfr(done, function(name) {
  v <- variants[[name]]
  class_dir <- file.path(out_root, ifelse(v$step == '03', 'default', name))
  evaluate_variant(file.path(out_root, name), ref_class, ref_dates,
                   class_dir = class_dir) %>%
    mutate(variant = name, group = v$group, parameter = v$parameter,
           value = v$value, default = v$default, step = v$step, .before = 1)
})
write_csv(per_track, file.path(res_dir, 'sens_per_track.csv'))


# =============================================================================
# 3. SUMMARY OF EACH VARIANT ----
# =============================================================================

base <- per_track %>% filter(variant == 'default') %>%
  select(ID, class_0 = class_agreement, within3_0 = prop_within3,
         mae_0 = mean_abs_error)

summary_tab <- per_track %>%
  left_join(base, by = 'ID') %>%
  group_by(variant, group, parameter, value, default, step) %>%
  summarise(
    # (new names: a summary with the same name as a column would replace it
    # in the next lines)
    n_tracks            = n(),
    class_agreement_pct = round(100 * weighted.mean(class_agreement, n_positions,
                                                    na.rm = TRUE), 1),
    class_median_track  = round(100 * median(class_agreement, na.rm = TRUE), 1),
    tracks_below_80     = sum(class_agreement < 0.8, na.rm = TRUE),
    dates_within3       = round(100 * mean(prop_within3, na.rm = TRUE), 1),
    dates_mean_abs_err  = round(mean(mean_abs_error, na.rm = TRUE), 1),
    dates_missing       = sum(n_missing, na.rm = TRUE),
    dates_check         = round(100 * mean(dates_check == 'Check', na.rm = TRUE), 1),
    segments_check      = round(100 * weighted.mean(prop_check, n_segments,
                                                    na.rm = TRUE), 1),
    n_segments_mean     = round(mean(n_segments, na.rm = TRUE), 1),
    run_time_median_s   = median(run_time_s, na.rm = TRUE),
    # Paired differences with the main analysis (percentage points / days)
    d_class      = 100 * mean(class_agreement - class_0, na.rm = TRUE),
    d_class_ci   = list(100 * boot_ci(class_agreement - class_0)),
    d_within3    = 100 * mean(prop_within3 - within3_0, na.rm = TRUE),
    d_within3_ci = list(100 * boot_ci(prop_within3 - within3_0)),
    d_mae        = mean(mean_abs_error - mae_0, na.rm = TRUE),
    d_mae_ci     = list(boot_ci(mean_abs_error - mae_0)),
    .groups = 'drop') %>%
  mutate(d_class_low    = map_dbl(d_class_ci, 1),   d_class_high    = map_dbl(d_class_ci, 2),
         d_within3_low  = map_dbl(d_within3_ci, 1), d_within3_high  = map_dbl(d_within3_ci, 2),
         d_mae_low      = map_dbl(d_mae_ci, 1),     d_mae_high      = map_dbl(d_mae_ci, 2)) %>%
  select(-ends_with('_ci')) %>%
  mutate(across(starts_with('d_'), ~ round(.x, 2)),
         variant = factor(variant, levels = names(variants))) %>%
  arrange(variant)

# Run time is not comparable when the segmentation was read from the cache
from_cache <- map_lgl(as.character(summary_tab$variant), function(name) {
  f <- file.path(out_root, name, 'Summary_segmentation.csv')
  variants[[name]]$step == '03' ||
    (file.exists(f) && any(read_csv(f, show_col_types = FALSE)$segclust_from_cache))
})
summary_tab$run_time_median_s[from_cache] <- NA

write_csv(summary_tab, file.path(res_dir, 'sens_summary.csv'))
print(summary_tab %>% select(variant, value, class_agreement_pct, dates_within3,
                             dates_mean_abs_err, segments_check, d_class,
                             d_within3), n = Inf)


# =============================================================================
# 4. EACH DATE SEPARATELY ----
# =============================================================================
# Bias (median difference) and proportion within 3 days for each date: some
# parameters (e.g. KDE level, departure rule) change only some dates.

dates_tab <- per_track %>%
  select(variant, parameter, value, ID, starts_with('err_')) %>%
  pivot_longer(starts_with('err_'), names_to = 'date', values_to = 'error',
               names_prefix = 'err_') %>%
  group_by(variant, parameter, value, date) %>%
  summarise(n = sum(!is.na(error)),
            bias_median = median(error, na.rm = TRUE),
            mean_abs_error = round(mean(abs(error), na.rm = TRUE), 1),
            within3 = round(100 * mean(abs(error) <= 3, na.rm = TRUE), 1),
            .groups = 'drop') %>%
  mutate(date = factor(date, levels = date_names))
write_csv(dates_tab, file.path(res_dir, 'sens_dates.csv'))


# =============================================================================
# 5. OPTIONAL: SECOND ANALYST ----
# =============================================================================
# Only the tracks with dates of both analysts. For each variant and each
# analyst: proportion of dates within 3 days. The row 'Analyst vs reference'
# is the variability between the two people.

if (!is.null(other)) {
  common <- intersect(other$ID, ref_dates$ID)
  cat('\nSecond analyst:', length(common), 'tracks in common\n')

  compare_dates <- function(a, b) {
    j <- inner_join(a, b, by = 'ID', suffix = c('_a', '_b'))
    e <- sapply(date_names, function(n) as.numeric(j[[paste0(n, '_a')]] -
                                                     j[[paste0(n, '_b')]]))
    e <- matrix(e, ncol = length(date_names), dimnames = list(NULL, date_names))
    tibble(n_tracks = nrow(j), n_dates = sum(!is.na(e)),
           within3 = round(100 * mean(abs(e) <= 3, na.rm = TRUE), 1),
           mean_abs_error = round(mean(abs(e), na.rm = TRUE), 1),
           bias_median = median(e, na.rm = TRUE))
  }

  other_tab <- map_dfr(unique(other$Analyst), function(an) {
    oa <- other %>% filter(Analyst == an, ID %in% common) %>% select(-Analyst)
    between_people <- compare_dates(oa, ref_dates %>% filter(ID %in% common)) %>%
      mutate(variant = 'Analyst vs reference', .before = 1)
    vs_variants <- map_dfr(done, function(name) {
      ph <- read_csv(file.path(out_root, name, 'Phenology.csv'),
                     show_col_types = FALSE, col_types = cols(ID = 'c'))
      auto <- main_dates(ph, suffix = '_KDE') %>% filter(ID %in% common)
      bind_rows(
        compare_dates(auto, ref_dates %>% filter(ID %in% common)) %>%
          mutate(compared_with = 'Reference'),
        compare_dates(auto, oa) %>% mutate(compared_with = 'Analyst')) %>%
        mutate(variant = name, .before = 1)
    })
    bind_rows(between_people, vs_variants) %>% mutate(Analyst = an, .before = 1)
  })
  write_csv(other_tab, file.path(res_dir, 'other_analyst.csv'))
  print(other_tab %>% filter(variant %in% c('Analyst vs reference', 'default')))
}
