# =============================================================================
# SENSITIVITY TESTS: LIST OF VARIANTS
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Each variant changes ONE parameter (or one pair of related
#              thresholds); all the other parameters keep the values of
#              parameters.R. 'default' is the main analysis, used as the
#              baseline of all comparisons.
#
#   group     'Requested' (asked by the reviewers) or 'Additional'
#   parameter name of the parameter in the tables and figures
#   value     value of the parameter in this variant (text)
#   default   value of the parameter in the main analysis (text)
#   step      '01' the variant changes the segmentation or the classification:
#                  scripts 01 and 03 are run again
#             '03' the variant only changes the phenological dates: only
#                  script 03 is run, on the classification of 'default'
#   changes   values that replace those of parameters.R (same names; the
#             values of 'adv' go inside adv = list(...)).
#             'positions_file' is the file with the positions (by default
#             data/All_positions.csv).
#
# Notes:
# - segclust2d does not accept segments shorter than 5 positions (2.5 days with
#   2 positions per day), so the minimum segment duration can only be tested
#   with longer values.
# - The minimum duration of a stage was a separate parameter (min_stage_days).
#   It had no effect on the dates (1.5 to 20 days), because a stage cannot be
#   shorter than a segment and W-mig-W sequences are already joined. It was
#   removed: the minimum stage length is now the minimum segment length.
# =============================================================================

variant <- function(group, parameter, value, default, step, changes = list()) {
  list(group = group, parameter = parameter, value = as.character(value),
       default = as.character(default), step = step, changes = changes)
}

# Positions with other treatments of the colony days (created with
# analyses/sensitivity/01_prepare_cdays_positions.R)
cdays_file <- function(name) {
  file.path('data', 'sensitivity', paste0('All_positions_cdays_', name, '.csv'))
}

variants <- list(

  # --- Main analysis -----------------------------------------------------------
  default = variant('Default', 'Default', '-', '-', '01'),

  # --- Requested by the reviewers ------------------------------------------------
  kmax_25      = variant('Requested', 'Kmax', 25, 15, '01', list(Kmax = 25)),
  kmax_35      = variant('Requested', 'Kmax', 35, 15, '01', list(Kmax = 35)),

  ncluster_3_6  = variant('Requested', 'ncluster', '3:6',  '5:10', '01',
                          list(ncluster = 3:6)),
  ncluster_5_15 = variant('Requested', 'ncluster', '5:15', '5:10', '01',
                          list(ncluster = 5:15)),

  lmin_3.5     = variant('Requested', 'Minimum segment duration (days)', 3.5, 2.5,
                         '01', list(min_segment_days = 3.5)),
  lmin_5       = variant('Requested', 'Minimum segment duration (days)', 5, 2.5,
                         '01', list(min_segment_days = 5)),

  maxstep_2220 = variant('Requested', 'Maximum step (km)', 2220, 4440, '01',
                         list(adv = list(max_step = 2220))),
  maxstep_none = variant('Requested', 'Maximum step (km)', 'None', 4440, '01',
                         list(adv = list(max_step = Inf))),

  eqwin_10     = variant('Requested', 'Equinox window (days)', 10, 20, '01',
                         list(eq_window_days = 10)),
  eqwin_30     = variant('Requested', 'Equinox window (days)', 30, 20, '01',
                         list(eq_window_days = 30)),

  cdays_keep   = variant('Requested', 'Colony days', 'Kept', 'Random (100 km)',
                         '01', list(positions_file = cdays_file('keep'))),
  cdays_remove = variant('Requested', 'Colony days', 'Removed', 'Random (100 km)',
                         '01', list(positions_file = cdays_file('remove'))),
  cdays_anchor = variant('Requested', 'Colony days', 'Anchored', 'Random (100 km)',
                         '01', list(positions_file = cdays_file('anchor'))),

  kde_75       = variant('Requested', 'KDE level (%)', 75, 90, '03',
                         list(kde_level = 75)),
  kde_95       = variant('Requested', 'KDE level (%)', 95, 90, '03',
                         list(kde_level = 95)),

  bw_0.5       = variant('Requested', 'KDE bandwidth (x href)', 0.5, 1, '03',
                         list(adv = list(kde_bandwidth_mult = 0.5))),
  bw_1.5       = variant('Requested', 'KDE bandwidth (x href)', 1.5, 1, '03',
                         list(adv = list(kde_bandwidth_mult = 1.5))),

  # --- Additional tests ------------------------------------------------------------
  # Thresholds of the agreement score (straightness and ang_vector together,
  # because they have the same threshold)
  thr_sa_0.4   = variant('Additional', 'Straightness and ang_vector threshold',
                         0.4, 0.5, '01',
                         list(adv = list(score_straightness = 0.4,
                                         score_ang_vector = 0.4))),
  thr_sa_0.6   = variant('Additional', 'Straightness and ang_vector threshold',
                         0.6, 0.5, '01',
                         list(adv = list(score_straightness = 0.6,
                                         score_ang_vector = 0.6))),
  thr_sin_0.8  = variant('Additional', 'Sinuosity threshold', 0.8, 1, '01',
                         list(adv = list(score_sinuosity = 0.8))),
  thr_sin_1.2  = variant('Additional', 'Sinuosity threshold', 1.2, 1, '01',
                         list(adv = list(score_sinuosity = 1.2))),
  r5_0.65      = variant('Additional', 'ang_vector of rule R5', 0.65, 0.75, '01',
                         list(adv = list(breeding_ang_vector = 0.65))),
  r5_0.85      = variant('Additional', 'ang_vector of rule R5', 0.85, 0.75, '01',
                         list(adv = list(breeding_ang_vector = 0.85))),

  cdays_r50    = variant('Additional', 'Colony days radius (km)', 50, 100, '01',
                         list(positions_file = cdays_file('random50'))),
  cdays_r200   = variant('Additional', 'Colony days radius (km)', 200, 100, '01',
                         list(positions_file = cdays_file('random200'))),
  # Other random points within 100 km (another seed of the random numbers)
  cdays_seed2  = variant('Additional', 'Colony days random seed', 2, 1, '01',
                         list(positions_file = cdays_file('seed2'))),
  cdays_seed3  = variant('Additional', 'Colony days random seed', 3, 1, '01',
                         list(positions_file = cdays_file('seed3'))),

  dep_first    = variant('Additional', 'Departure', 'First exit', 'Last exit',
                         '03', list(adv = list(departure_rule = 'first_exit')))

)

# Treatments of the colony days used above (for 01_prepare_cdays_positions.R)
cdays_variants <- list(
  keep      = list(mode = 'keep',   radius = 100),
  remove    = list(mode = 'remove', radius = 100),
  anchor    = list(mode = 'anchor', radius = 100),
  random50  = list(mode = 'random', radius = 50),
  random200 = list(mode = 'random', radius = 200),
  seed2     = list(mode = 'random', radius = 100, seed = 2),
  seed3     = list(mode = 'random', radius = 100, seed = 3)
)
