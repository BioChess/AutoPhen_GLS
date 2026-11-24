# =============================================================================
# BASELINE METHODS 02: COMPARISON WITH THE WORKFLOW
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Compares the three baselines and the workflow (AutoPhen_GLS,
#              fully automated output) with the truth of the simulations and
#              with the supervised reference of the real tracks, with the same
#              measures:
#     class_agreement  proportion of positions with the true (reference) class
#                      (breeding / migration / wintering)
#     err_<date>       date error (days; estimated - true)
#     within3          proportion of the four dates within 3 days; a date
#                      that the method could not estimate counts as wrong
#   Summary (per dataset, observation error and method): mean class agreement,
#   dates within 1, 3 and 7 days, median absolute error, missing dates, failed
#   tracks, and the mean paired difference with the workflow (method -
#   workflow, the same tracks) with a bootstrap 95% confidence interval.
#
#   As in the evaluation of the simulations (simulations/03), the true class
#   of a movement between two wintering stages is wintering.
#
# Input:  output/baselines/<dataset>/ (01_run_baselines.R)
#         output/simulations/<gls|none>/ and output/sensitivity/default/
#         (results of the workflow)
#         simulations/data/ (truth); data/reference/ (supervised reference)
# Output: output/baselines/_results/baselines_per_track.csv
#         output/baselines/_results/baselines_summary.csv
#         output/baselines/_results/baselines_by_distance.csv (simulations)
#         output/baselines/_results/fig_baselines.png / .pdf
# =============================================================================

pacman::p_load('tidyverse')
source('analyses/sensitivity/sensitivity_functions.R')   # three_classes(), main_dates(), boot_ci()
source('analyses/baselines/baseline_functions.R')        # final_classes(), dates_from_classes()
env <- new.env()
sys.source(file.path('functions', 'phenology_dates_functions.R'), envir = env)
set.seed(1)

res_dir <- file.path('output', 'baselines', '_results')
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
methods <- c('segclust_simple', 'migrateR', 'hmm')
method_labels <- c(autophen = 'AutoPhen_GLS (this workflow)',
                   autophen_transition = 'AutoPhen_GLS classes, dates from class changes',
                   segclust_simple = 'segclust2d + simple rule',
                   migrateR = 'migrateR (NSD models)',
                   hmm = 'Hidden Markov model (3 states)')


# =============================================================================
# 1. TRUTH AND REFERENCE ----
# =============================================================================

# Movement between two wintering stages = wintering (as in the workflow)
join_W <- function(d) {
  d %>% arrange(ID, date) %>% group_by(ID) %>%
    mutate(prev = if_else(ref == 'mig', NA_character_, ref), nxt = prev) %>%
    fill(prev, .direction = 'down') %>% fill(nxt, .direction = 'up') %>%
    mutate(ref = if_else(ref == 'mig' & prev %in% 'W' & nxt %in% 'W', 'W', ref)) %>%
    ungroup() %>% select(ID, date, ref)
}

refs <- list()
if (file.exists(file.path('simulations', 'data', 'sim_truth_positions.csv'))) {
  sim_md <- read_csv(file.path('simulations', 'data', 'sim_metadata.csv'),
                     show_col_types = FALSE, col_types = cols(geo_trip = 'c')) %>%
    select(ID = geo_trip, distance, movement, error)
  tp <- read_csv(file.path('simulations', 'data', 'sim_truth_positions.csv'),
                 show_col_types = FALSE, col_types = cols(geo_trip = 'c'),
                 col_select = c(geo_trip, Date_Time, stage)) %>%
    transmute(ID = geo_trip, date = Date_Time, ref = stage) %>% join_W()
  td <- read_csv(file.path('simulations', 'data', 'sim_truth_dates.csv'),
                 show_col_types = FALSE) %>%
    transmute(ID = sim_id, B1_dep, W_arr, W_dep, B2_arr)
  refs$sim_gls  <- list(class = tp, dates = td,
                        ids = sim_md$ID[sim_md$error %in% c('gls', 'gls_eq')],
                        workflow = file.path('output', 'simulations', 'gls'))
  refs$sim_none <- list(class = tp, dates = td,
                        ids = sim_md$ID[sim_md$error == 'none'],
                        workflow = file.path('output', 'simulations', 'none'))
}
if (file.exists(file.path('data', 'reference', 'CAL_tracks_supervised.csv'))) {
  rc <- read_csv(file.path('data', 'reference', 'CAL_tracks_supervised.csv'),
                 show_col_types = FALSE,
                 col_types = cols(ID = 'c', date = col_datetime())) %>%
    transmute(ID, date, ref = three_classes(Breed2))
  rd <- read_csv(file.path('data', 'reference', 'CAL_phen_supervised.csv'),
                 show_col_types = FALSE, col_types = cols(ID = 'c')) %>% main_dates()
  refs$real <- list(class = rc, dates = rd, ids = unique(rc$ID),
                    workflow = file.path('output', 'sensitivity', 'default'))
}


# =============================================================================
# 2. PER-TRACK RESULTS ----
# =============================================================================

evaluate_method <- function(classes, dates, ref, ids) {
  cls <- inner_join(ref$class %>% filter(ID %in% ids), classes,
                    by = c('ID', 'date')) %>%
    group_by(ID) %>% summarise(class_agreement = mean(ref == class), .groups = 'drop')
  tru <- ref$dates %>% filter(ID %in% ids) %>% mutate(across(-ID, as.Date))
  est <- dates %>% mutate(ID = as.character(ID), across(-ID, as.Date))
  dts <- left_join(tru, est, by = 'ID', suffix = c('_true', '_est'))
  for (n in date_names) {
    dts[[paste0('err_', n)]] <- as.numeric(dts[[paste0(n, '_est')]] -
                                             dts[[paste0(n, '_true')]])
  }
  dts %>% select(ID, starts_with('err_')) %>%
    left_join(cls, by = 'ID') %>%
    mutate(done = ID %in% unique(classes$ID))
}

per_track <- list()
for (d in names(refs)) {
  ref <- refs[[d]]
  # Workflow
  wf <- ref$workflow
  if (file.exists(file.path(wf, 'Phenology.csv'))) {
    cl <- read_csv(file.path(wf, 'Phenology_tracks.csv'), show_col_types = FALSE,
                   col_types = cols(ID = 'c')) %>%
      transmute(ID, date, class = three_classes(Breed2))
    ph <- read_csv(file.path(wf, 'Phenology.csv'), show_col_types = FALSE,
                   col_types = cols(ID = 'c'))
    per_track[[paste(d, 'autophen')]] <- evaluate_method(
      cl, main_dates(ph, suffix = '_KDE'), ref, ref$ids) %>%
      mutate(dataset = d, method = 'autophen')
    # The same classes, with the dates of the baselines (first day of each
    # change of class) instead of the KDE: separates the effect of the
    # classification from the effect of the definition of the dates
    lmin <- readRDS(file.path(wf, 'parameters.rds'))$lmin
    cl_t <- cl %>% group_by(ID) %>% arrange(date, .by_group = TRUE) %>%
      mutate(class = final_classes(class, lmin, env)) %>% ungroup()
    dt_t <- cl_t %>% group_by(ID) %>%
      group_modify(~ as.data.frame(dates_from_classes(.x$date, .x$class))) %>% ungroup()
    per_track[[paste(d, 'autophen_transition')]] <- evaluate_method(
      cl_t, dt_t, ref, ref$ids) %>%
      mutate(dataset = d, method = 'autophen_transition')
  } else {
    cat('Results of the workflow not found for', d, '(', wf, ')\n')
  }
  # Baselines
  for (m in methods) {
    f <- file.path('output', 'baselines', d, paste0(m, '_dates.csv'))
    if (!file.exists(f)) { cat('Not found:', f, '\n'); next }
    cl <- read_csv(file.path('output', 'baselines', d, paste0(m, '_classes.csv')),
                   show_col_types = FALSE, col_types = cols(ID = 'c'))
    dt <- read_csv(f, show_col_types = FALSE, col_types = cols(ID = 'c'))
    per_track[[paste(d, m)]] <- evaluate_method(cl, dt, ref, ref$ids) %>%
      mutate(dataset = d, method = m)
  }
}
per_track <- bind_rows(per_track)
if (exists('sim_md')) per_track <- per_track %>% left_join(sim_md, by = 'ID')
if (!'error' %in% names(per_track)) per_track$error <- NA_character_
per_track <- per_track %>%
  mutate(error = ifelse(dataset == 'real', 'real (GLS)', error),
         within3 = rowMeans(across(starts_with('err_'),
                                   ~ replace_na(abs(.x) <= 3, FALSE))))
write_csv(per_track, file.path(res_dir, 'baselines_per_track.csv'))


# =============================================================================
# 3. SUMMARY ----
# =============================================================================

long_dates <- per_track %>%
  pivot_longer(starts_with('err_'), names_to = 'date', values_to = 'e')
summ <- function(g) {
  dates_s <- long_dates %>% group_by(across(all_of(g))) %>%
    summarise(n_dates = n(),
              missing_pct = round(100 * mean(is.na(e)), 1),
              within1_pct = round(100 * mean(replace_na(abs(e) <= 1, FALSE)), 1),
              within3_pct = round(100 * mean(replace_na(abs(e) <= 3, FALSE)), 1),
              within7_pct = round(100 * mean(replace_na(abs(e) <= 7, FALSE)), 1),
              abs_median  = median(abs(e), na.rm = TRUE), .groups = 'drop')
  cls_s <- per_track %>% group_by(across(all_of(g))) %>%
    summarise(n_tracks = n(), failed = sum(!done),
              class_agreement_pct = round(100 * mean(class_agreement, na.rm = TRUE), 1),
              .groups = 'drop')
  left_join(cls_s, dates_s, by = g)
}

# Paired differences with the workflow (same tracks)
paired <- function(g) {
  wf <- per_track %>% filter(method == 'autophen') %>%
    select(ID, dataset, ca_wf = class_agreement, w3_wf = within3)
  per_track %>% filter(!method %in% c('autophen')) %>%
    inner_join(wf, by = c('ID', 'dataset')) %>%
    group_by(across(all_of(g))) %>%
    summarise(
      class_diff = round(100 * mean(class_agreement - ca_wf, na.rm = TRUE), 1),
      class_diff_ci = paste(round(100 * boot_ci(class_agreement - ca_wf), 1), collapse = ' to '),
      within3_diff = round(100 * mean(within3 - w3_wf, na.rm = TRUE), 1),
      within3_diff_ci = paste(round(100 * boot_ci(within3 - w3_wf), 1), collapse = ' to '),
      .groups = 'drop')
}

g1 <- c('dataset', 'error', 'method')
tab <- summ(g1) %>% left_join(paired(g1), by = g1) %>%
  mutate(method = factor(method, levels = names(method_labels))) %>%
  arrange(dataset, error, method)
write_csv(tab, file.path(res_dir, 'baselines_summary.csv'))
print(tab %>% select(dataset, error, method, failed, class_agreement_pct,
                     within3_pct, missing_pct, class_diff, within3_diff), n = Inf)

if (any(!is.na(per_track$distance))) {
  g2 <- c('error', 'distance', 'method')
  tab2 <- summ(g2) %>% filter(!is.na(distance)) %>%
    mutate(method = factor(method, levels = names(method_labels)),
           distance = factor(distance, levels = c('short', 'medium', 'long'))) %>%
    arrange(error, distance, method)
  write_csv(tab2, file.path(res_dir, 'baselines_by_distance.csv'))
}


# =============================================================================
# 4. FIGURE ----
# =============================================================================
# Mean (and bootstrap 95% CI) of the class agreement and of the proportion of
# dates within 3 days, per method and type of data

error_labels <- c(none = 'Simulated, no error', gls = 'Simulated, GLS-like error',
                  gls_eq = 'Simulated, GLS-like error,\nmigration at equinoxes',
                  `real (GLS)` = 'Calonectris GLS tracks\n(supervised reference)')
fig_dat <- per_track %>%
  pivot_longer(c(class_agreement, within3), names_to = 'measure', values_to = 'v') %>%
  group_by(error, method, measure) %>%
  summarise(mean = mean(v, na.rm = TRUE),
            lo = boot_ci(v)[1], hi = boot_ci(v)[2], .groups = 'drop') %>%
  mutate(method = factor(method_labels[method], levels = rev(method_labels)),
         error = factor(error_labels[error], levels = error_labels),
         measure = recode(measure, class_agreement = 'Positions correctly classified',
                          within3 = 'Dates within 3 days'))

fig <- ggplot(fig_dat, aes(x = 100 * mean, y = method, colour = method == method_labels[['autophen']])) +
  geom_errorbarh(aes(xmin = 100 * lo, xmax = 100 * hi), height = 0.25) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = 'black', `FALSE` = 'grey55'), guide = 'none') +
  facet_grid(error ~ measure) +
  labs(x = 'Percentage (mean per track and 95% CI)', y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'),
        strip.text.y = element_text(angle = 0, hjust = 0), panel.grid.minor = element_blank())
ggsave(file.path(res_dir, 'fig_baselines.png'), fig, width = 22, height = 16,
       units = 'cm', dpi = 300)
ggsave(file.path(res_dir, 'fig_baselines.pdf'), fig, width = 22, height = 16, units = 'cm')
cat('Results in', res_dir, '\n')
