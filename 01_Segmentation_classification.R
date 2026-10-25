# =============================================================================
# 01. SEGMENTATION, CLUSTERING AND CLASSIFICATION OF TRACKS (all tracks)
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Runs step 1 of the workflow for every track in the metadata
#              file: segmentation and clustering of the track (segclust2d),
#              classification of segments into Breeding / Migration /
#              Wintering, and classification of each position.
#              The same steps are explained one by one, for one track, in
#              01_Segmentation_classification.qmd.
#
# Input:  positions (one row per position: geo_trip, Date_Time, lon, lat)
#         metadata  (one row per track: geo_trip, Ring, Sp, Colony,
#                    Colony.lon, Colony.lat)
# Output (in output/, or output/full/ for the full dataset):
#         Phenology_tracks.csv             (one row per position)
#         output/Segmentation_results.csv  (one row per segment)
#         output/Summary_segmentation.csv  (one row per track)
#         output/plots/segmentation/*.png  (optional maps)
# =============================================================================


# =============================================================================
# 1. PACKAGES, FUNCTIONS AND PARAMETERS ----
# =============================================================================

pacman::p_load('tidyverse', 'lubridate', 'readr', 'sf',
               'adehabitatLT', 'adehabitatHR', 'segclust2d', 'trajr',
               'recipes', 'furrr', 'units', 'ggplot2', 'patchwork',
               'rnaturalearth')

source('functions/segmentation_classification_functions.R')
source('parameters.R')


# =============================================================================
# 2. SETTINGS OF THIS SCRIPT ----
# =============================================================================

# Data to analyse:
#   'GLS'  example data: Cory's shearwater GLS tracks (input/)
#   'GPS'  example data: Audouin's gull GPS track (input/)
#   'full' full dataset of the paper (data/, not shared; created with
#          analyses/prepare_full_dataset.R). Results are saved in output/full/
dataset <- 'GLS'
if (dataset == 'full') {
  positions_file <- 'data/All_positions.csv'
  metadata_file  <- 'data/All_metadata.csv'
  out_dir        <- 'output/full'
} else {
  positions_file <- paste0('input/', dataset, '_positions_example.csv')
  metadata_file  <- paste0('input/', dataset, '_metadata_example.csv')
  out_dir        <- 'output'
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save a map of the segmentation and classification of each track?
save_plots <- TRUE
plots_dir  <- file.path(out_dir, 'plots', 'segmentation')


# =============================================================================
# 3. READ DATA ----
# =============================================================================

metadata <- read_csv(metadata_file, show_col_types = FALSE,
                     col_types = cols(geo_trip = 'c'))
tracks   <- read_csv(positions_file, show_col_types = FALSE,
                     col_types = cols(geo_trip = 'c'))

if (save_plots) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 4. LOOP OVER TRACKS ----
# =============================================================================
# Each track is processed inside tryCatch(): if one track fails, the error is
# reported and the loop continues with the next track.

seg_list <- list(); summary_list <- list(); track_list <- list()
failed   <- data.frame(ID = character(), error = character())

for (i in seq_len(nrow(metadata))) {

  info  <- metadata[i, ]
  track <- tracks %>% filter(geo_trip == info$geo_trip)
  cat('\n[', i, '/', nrow(metadata), '] Track', info$geo_trip, '-',
      nrow(track), 'positions\n')
  t0 <- Sys.time()

  res <- tryCatch(
    classify_track(track = track, info = info,
                   Kmax = Kmax, lmin = lmin, ncluster = ncluster,
                   eq_window_days = eq_window_days, adv = adv,
                   make_plots = save_plots),
    error = function(e) {
      message('  ERROR in track ', info$geo_trip, ': ', conditionMessage(e))
      failed <<- rbind(failed, data.frame(ID = as.character(info$geo_trip),
                                          error = conditionMessage(e)))
      NULL
    })
  if (is.null(res)) next

  id <- as.character(info$geo_trip)
  seg_list[[id]]     <- res$segments
  summary_list[[id]] <- res$summary %>%
    mutate(run_time_s = round(as.numeric(difftime(Sys.time(), t0,
                                                  units = 'secs'))))
  track_list[[id]]   <- res$track

  if (save_plots && !is.null(res$plot)) {
    ggsave(file.path(plots_dir, paste0(id, '_segmentation.png')), res$plot,
           width = 30, height = 15, units = 'cm', dpi = 150)
  }
  cat('  Done in', summary_list[[id]]$run_time_s, 's\n')
}


# =============================================================================
# 5. SAVE RESULTS ----
# =============================================================================

seg_all     <- bind_rows(lapply(seg_list, function(d) mutate(d, ID = as.character(ID))))
summary_all <- bind_rows(lapply(summary_list, function(d) mutate(d, ID = as.character(ID))))
track_all   <- bind_rows(lapply(track_list, function(d) mutate(d, ID = as.character(ID))))

write_csv(seg_all,     file.path(out_dir, 'Segmentation_results.csv'))
write_csv(summary_all, file.path(out_dir, 'Summary_segmentation.csv'))
write_csv(track_all,   file.path(out_dir, 'Phenology_tracks.csv'))
if (nrow(failed) > 0) write_csv(failed, file.path(out_dir, 'Failed_tracks_01.csv'))

cat('\nFinished:', length(track_list), 'tracks classified,', nrow(failed),
    'failed.\n')
