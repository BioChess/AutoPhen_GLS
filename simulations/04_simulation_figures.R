# =============================================================================
# SIMULATIONS 04: FIGURES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Figures of the simulation study.
#   fig_sim_examples_<distance>  one example track per movement type (columns)
#                                and observation error (rows)
#   fig_sim_parameters           values of the parameters of the simulated
#                                tracks (km and days), by distance class
#   fig_sim_date_errors          error of each date (estimated - true, days)
#                                by scenario
#   fig_sim_ratios               results of each track as a function of the
#                                ratio of migration distance to breeding-season
#                                range and to observation error
#
# Input:  simulations/data/ (01_simulate_tracks.R)
#         output/simulations/_results/sim_per_track.csv (03_evaluate_simulations.R)
# Output: output/simulations/_results/fig_sim_*.png and .pdf
# =============================================================================

pacman::p_load('tidyverse', 'patchwork', 'scales')
source('simulations/simulation_functions.R')   # labels, colours, plot_sim_tracks()

sim_dir <- file.path('simulations', 'data')
res_dir <- file.path('output', 'simulations', '_results')
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

metadata    <- read_csv(file.path(sim_dir, 'sim_metadata.csv'), show_col_types = FALSE)
truth_dates <- read_csv(file.path(sim_dir, 'sim_truth_dates.csv'), show_col_types = FALSE)
truth_pos   <- read_csv(file.path(sim_dir, 'sim_truth_positions.csv'), show_col_types = FALSE)
p           <- readRDS(file.path(sim_dir, 'sim_settings.rds'))

save_fig <- function(g, name, width, height) {
  ggsave(file.path(res_dir, paste0(name, '.png')), g, width = width, height = height,
         units = 'cm', dpi = 300)
  ggsave(file.path(res_dir, paste0(name, '.pdf')), g, width = width, height = height,
         units = 'cm')
}
distance_fills <- setNames(c('#c6dbef', '#6baed6', '#2171b5'), distance_labels)


# =============================================================================
# 1. EXAMPLE TRACKS ----
# =============================================================================
# First track of each scenario (rep = 1)

first_ids <- truth_dates %>%
  group_by(distance, movement, error) %>% slice(1) %>% ungroup()
for (d in c('short', 'medium', 'long')) {
  ids <- first_ids %>% filter(distance == d) %>% pull(sim_id)
  dat <- truth_pos %>% filter(geo_trip %in% ids) %>%
    left_join(metadata %>% select(geo_trip, movement, error), by = 'geo_trip')
  save_fig(plot_sim_tracks(dat, title = distance_labels[[d]]),
           paste0('fig_sim_examples_', d), 24, 19)
}


# =============================================================================
# 2. PARAMETERS OF THE SIMULATED TRACKS ----
# =============================================================================
# Real scales of the simulations: what each distance class means in km and
# days, and the ratios of migration distance to breeding-season range and to
# observation error

eq_dist <- function(d, md) {
  abs(as.numeric(difftime(d, as.POSIXct(paste0(format(d, '%Y'), md), tz = 'UTC'),
                          units = 'days')))
}
per_track_sim <- truth_pos %>%
  mutate(r_true = sqrt(x_true_km^2 + y_true_km^2),
         near_eq = eq_dist(Date_Time, '-09-22') <= p$error$gls$eq_days |
                   eq_dist(Date_Time, '-03-20') <= p$error$gls$eq_days) %>%
  group_by(sim_id = geo_trip) %>%
  summarise(breeding_range_km = unname(quantile(r_true[stage == 'B'], 0.95)),
            obs_error_km = mean(error_km[!near_eq]),
            .groups = 'drop')

pars <- truth_dates %>%
  left_join(per_track_sim, by = 'sim_id') %>%
  mutate(out_days    = as.numeric(W_arr - B1_dep),
         return_days = as.numeric(B2_arr - W_dep),
         winter_days = as.numeric(W_dep - W_arr),
         tau_v_hours = tau_v_days * 24,
         ratio_range = dist_km / breeding_range_km,
         ratio_error = ifelse(error == 'none', NA, dist_km / obs_error_km))
write_csv(pars, file.path(res_dir, 'sim_track_parameters.csv'))

par_labels <- c(
  dist_km = 'Distance to wintering area (km)',
  speed_out = 'Migration speed, outward (km/day)',
  range_radius_km = 'Radius of the stationary areas (km, 95%)',
  breeding_range_km = 'Breeding-season range (km, 95%)',
  tau_z_days = 'Position autocorrelation, tau.z (days)',
  tau_v_hours = 'Velocity autocorrelation, tau.v (hours)',
  out_days = 'Outward migration (days)',
  return_days = 'Return migration (days)',
  winter_days = 'Wintering period (days)',
  ratio_range = 'Distance / breeding-season range',
  ratio_error = 'Distance / observation error (GLS)')

pars_long <- pars %>%
  select(sim_id, distance, movement, all_of(names(par_labels))) %>%
  pivot_longer(all_of(names(par_labels)), names_to = 'parameter') %>%
  filter(!is.na(value)) %>%
  mutate(parameter = factor(par_labels[parameter], levels = par_labels)) %>%
  label_scenarios()

g_par <- ggplot(pars_long, aes(distance, value, fill = distance)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.4) +
  scale_fill_manual(values = distance_fills, name = 'Migration distance') +
  scale_x_discrete(labels = c('Short', 'Medium', 'Long')) +
  facet_wrap(~ parameter, scales = 'free_y', ncol = 4,
             labeller = label_wrap_gen(28)) +
  labs(x = NULL, y = NULL) +
  theme_sim
save_fig(g_par, 'fig_sim_parameters', 22, 18)


# =============================================================================
# 3. RESULTS ----
# =============================================================================

per_track <- read_csv(file.path(res_dir, 'sim_per_track.csv'), show_col_types = FALSE) %>%
  left_join(pars %>% select(ID = sim_id, ratio_range, ratio_error, breeding_range_km),
            by = 'ID')
date_cols <- c(err_B1_dep = 'Departure from\nbreeding area', err_W_arr = 'Arrival at\nwintering area',
               err_W_dep = 'Departure from\nwintering area', err_B2_arr = 'Arrival at\nbreeding area')

# --- 3.1. Error of each date by scenario ------------------------------------------
errs <- per_track %>%
  select(ID, distance, movement, error, all_of(names(date_cols))) %>%
  pivot_longer(all_of(names(date_cols)), names_to = 'date', values_to = 'diff') %>%
  filter(!is.na(diff)) %>%
  mutate(date = factor(date_cols[date], levels = date_cols)) %>%
  label_scenarios()

g_err <- ggplot(errs, aes(movement, diff, fill = distance)) +
  geom_hline(yintercept = 0, colour = 'grey30', linewidth = 0.3) +
  geom_hline(yintercept = c(-3, 3), colour = 'grey60', linewidth = 0.3, linetype = 2) +
  geom_boxplot(outlier.size = 0.6, outlier.shape = 21, linewidth = 0.3,
               position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = distance_fills, name = 'Migration distance') +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(-100, -30, -10, -3, 0, 3, 10, 30, 100)) +
  scale_x_discrete(labels = label_wrap(10)) +
  facet_grid(error ~ date) +
  labs(x = NULL, y = 'Estimated - true date (days)',
       caption = 'Dashed lines: ±3 days. Pseudo-logarithmic scale.') +
  theme_sim
save_fig(g_err, 'fig_sim_date_errors', 24, 18)

# --- 3.2. Results as a function of the ratios ---------------------------------------
res <- per_track %>%
  mutate(date_error = rowMeans(abs(across(all_of(names(date_cols)))), na.rm = TRUE),
         date_error = ifelse(is.nan(date_error), NA, date_error)) %>%
  label_scenarios()
movement_colours <- setNames(c('grey30', 'orange2', '#7b3294', 'steelblue'), movement_labels)

ratio_plot <- function(d, x, xlab, y, ylab, ytrans = 'identity', ybreaks = waiver()) {
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_point(aes(fill = movement), shape = 21, colour = 'black', stroke = 0.15,
               size = 1.6, alpha = 0.8) +
    geom_smooth(method = 'loess', formula = y ~ x, se = TRUE, colour = 'black',
                linewidth = 0.5, fill = 'grey80') +
    scale_x_log10() +
    scale_y_continuous(trans = ytrans, breaks = ybreaks) +
    scale_fill_manual(values = movement_colours, name = 'Movement') +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    facet_wrap(~ error, nrow = 1) +
    labs(x = xlab, y = ylab) +
    theme_sim
}
g_ratio <- (ratio_plot(res, 'ratio_range', NULL, 'date_error',
                       'Mean absolute date error\n(days per track)',
                       ytrans = pseudo_log_trans(sigma = 1, base = 10),
                       ybreaks = c(0, 1, 3, 10, 30, 100)) /
            ratio_plot(res, 'ratio_range',
                       'Migration distance / breeding-season range (log scale)',
                       'class_agreement', 'Positions correctly classified\n(proportion per track)')) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
save_fig(g_ratio, 'fig_sim_ratios', 24, 18)

g_ratio_err <- ratio_plot(res %>% filter(!is.na(ratio_error)), 'ratio_error',
                          'Migration distance / mean observation error (log scale)',
                          'date_error', 'Mean absolute date error\n(days per track)',
                          ytrans = pseudo_log_trans(sigma = 1, base = 10),
                          ybreaks = c(0, 1, 3, 10, 30, 100))
save_fig(g_ratio_err, 'fig_sim_ratio_error', 18, 10)
