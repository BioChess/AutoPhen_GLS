# =============================================================================
# PARAMETERS OF THE AutoPhen_GLS WORKFLOW
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: All settings used by scripts 01-03 are defined here.
#
#   PART 1 - USER PARAMETERS
#     A few settings that users may need to adapt to their data or species.
#     Default values were chosen for seabirds tracked with GLS (two positions
#     per day). Please read the comment of each one before changing it.
#
#   PART 2 - ADVANCED PARAMETERS
#     Thresholds and rules that are part of the classification method.
#     They are listed here to make the method transparent, but they are NOT
#     meant to be changed by most users. Changing them changes the method.
#
# IDs in brackets (e.g. [P05]) refer to the parameter inventory (Table SX).
#
# How to use: source('parameters.R') at the start of scripts 01-03.
# =============================================================================


# =============================================================================
# PART 1 - USER PARAMETERS ----
# =============================================================================

# --- 1.1. Type of data --------------------------------------------------------
# "GLS": light-level geolocator positions. Rules to deal with the low accuracy
#        of latitude around the equinoxes are applied.
# "GPS": GPS (or other accurate) positions. Equinox rules are not applied.
data_type <- "GLS"

# --- 1.2. Positions per day ---------------------------------------------------
# Number of positions per day in the input data (GLS: 2, one per twilight).
# Used to convert durations in days into numbers of positions.
# For GPS data, we recommend resampling to 2 positions per day
# (see analyses/prepare_GPS_example.R).
positions_per_day <- 2

# --- 1.3. Equinox window (only for GLS data) [P05] -----------------------------
# Number of days before and after each equinox (around 21 March and
# 21 September) when GLS latitude is unreliable. Segments in this window are
# flagged and may be reclassified.
eq_window_days <- 20

# --- 1.4. Segmentation [P09, P10, P11] -----------------------------------------
# Kmax: maximum number of segments per track.
#   Increase it for long tracks or tracks with many movement phases.
Kmax <- 25

# ncluster: range of numbers of clusters tested by segclust2d. The best value
#   is selected by BIC.
ncluster <- 5:10

# min_segment_days: minimum duration of a segment (days). It is converted to
#   positions (lmin) with positions_per_day. segclust2d needs at least 5
#   positions per segment, so lmin is never lower than 5.
min_segment_days <- 2.5

# --- 1.5. Phenological stages [P36] -------------------------------------------
# min_stage_days: stages shorter than this (days) that are between two equal
#   stages are joined to them (e.g. a 2-day 'migration' inside the wintering
#   period becomes wintering).
min_stage_days <- 2.5

# --- 1.6. Breeding and wintering areas [P41] ----------------------------------
# kde_level: percentage of the kernel density estimation (KDE) that defines the
#   breeding and wintering areas. Arrival = first position inside the KDE;
#   departure = first position outside the KDE after the last position inside.
kde_level <- 90


# =============================================================================
# PART 2 - ADVANCED PARAMETERS (part of the method; do not change) ----
# =============================================================================
# Distances are in km: positions are projected to an azimuthal equidistant
# projection centred on the colony [P02].

adv <- list(

  # --- 2.1. Preparation of the track ------------------------------------------
  max_step          = 4440, # [P03] (km) positions with a longer step to the
                            #   next position are removed (outliers).
                            #   Equivalent to the 40 degrees used before.
                            #   TO REVIEW: a maximum speed may be better.
  max_repeated      = 4,    # [P07] more repeated latitudes than this → removed
  edge_positions    = 5,    # [P08] positions checked for repetitions at the
                            #   start and end of the track

  # --- 2.2. Distance threshold breeding / non-breeding [P19] ------------------
  nscd_Kmax         = 3,    # segments for the segmentation of NSCD
  nscd_ncluster     = 3,    # clusters for the segmentation of NSCD

  # --- 2.3. Trajectory metrics [P15, P16] -------------------------------------
  core_metrics      = c("straightness", "sinuosity", "ang_vector", "mu.dist"),
  extra_metrics     = c("sd.dist", "Emax", "sq_disp",
                        "R2n_sd", "dc_mean", "dc_sd"),
  paired_metrics    = list(mu_xy  = c("mu.x", "mu.y"),
                           sd_xy  = c("sd.x", "sd.y"),
                           mu_dxy = c("mu.dx", "mu.dy"),
                           sd_dxy = c("sd.dx", "sd.dy")),
  max_extra_metrics = 3,    # combinations of 1 to 3 extra metric groups

  # --- 2.4. Clustering of segments [P17] --------------------------------------
  clust_method      = "complete",  # hierarchical clustering linkage
  clust_k           = 2,           # clusters per clustering

  # --- 2.5. Classification rules [P20, P21] -----------------------------------
  short_straightness = 0.5,  # short stationary non-breeding segment with
                             #   higher straightness → migration
  breeding_ang_vector = 0.75, # migratory segment near the colony with lower
                             #   ang_vector → breeding

  # --- 2.6. Agreement score [P22, P24] ----------------------------------------
  # Agreement between the class of a segment and its movement metrics:
  # Breeding/wintering: straightness < 0.5, ang_vector < 0.5, sinuosity > 1,
  #   mean step length (mu.dist) <= 1.5 x median step length of the track
  # Migration: the opposite. Agreement = proportion of criteria met (0, 0.25,
  # 0.5, 0.75 or 1); 0 for impossible classes (see agreement_score()).
  # Sinuosity is dimensionless (see trajmetrics()): sinuosity > 1 means that
  # the mean cosine of turning angles is lower than 0.6.
  score_straightness = 0.5,
  score_ang_vector   = 0.5,
  score_sinuosity    = 1,
  score_speed        = 1.5,
  impossible_straightness = 0.2,  # migration with lower straightness and
  impossible_ang_vector   = 0.9,  #   lower ang_vector → agreement 0

  # --- 2.7. Equinox (only GLS) [P26, P27] -------------------------------------
  eq_sd_ratio_mig   = 2,    # migration → B/W if SD(y) > 2 · SD(x)
  eq_sd_limit       = 111,  # (km) breeding → migration if SD(x) or SD(y)
  eq_sd_ratio_breed = 4,    #   > 111 km and SD(x) > 4 · SD(y)

  # --- 2.8. Corrections [P28, P29, P31] ---------------------------------------
  corr_straightness = 0.2,
  corr_score        = 0.5,  # agreement (0-1)
  corr_sinuosity    = 1,
  # Fast wintering segments at the start or end of the wintering stage →
  # migration: mean step length > fast_w_ratio x that of the main wintering
  # segment and > fast_w_min (km between consecutive positions)
  fast_w_ratio      = 2.75,
  fast_w_min        = 250,

  # --- 2.9. Segments to review [P34] ------------------------------------------
  check_score       = 0.5,  # segments with lower agreement (0 or 0.25) are
                            #   flagged 'Check'
  w_far_km          = 1000, # wintering segments farther than this from the
                            #   main wintering segment are flagged 'Check'
                            #   (possible stopovers or migration loops)

  # --- 2.10. Kernel density estimation [P38, P40, P42] ------------------------
  kde_min_positions = 5,    # minimum positions to calculate a KDE
  kde_bandwidth     = "href",
  kde_fallback_grid = 1000,
  kde_fallback_extent = 4
)


# =============================================================================
# DERIVED VALUES (calculated from the parameters above; do not edit) ----
# =============================================================================

use_equinox <- data_type == "GLS"

# Minimum segment length in positions (segclust2d requires lmin >= 5)
lmin <- max(5, ceiling(min_segment_days * positions_per_day))

# Minimum stage length in positions
min_stage_length <- max(1, ceiling(min_stage_days * positions_per_day))
