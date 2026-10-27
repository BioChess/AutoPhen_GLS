# =============================================================================
# 03. PHENOLOGICAL DATES (all tracks)
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Runs step 3 of the workflow for every track classified in
#              step 1 (or corrected in step 2): delimits breeding and
#              wintering areas with kernel density estimation (KDE) and
#              extracts arrival and departure dates.
#              The same steps are explained one by one, for one track, in
#              03_Phenology_dates.qmd.
#
#              Definitions:
#                arrival   = first position inside the KDE
#                departure = first position outside the KDE after the last
#                            position inside it
#
# Input:  output/Phenology_tracks.csv (or _corrected.csv) and the metadata
# Output: output/Phenology.csv                 (one row per track)
#         output/plots/phenology/*.png         (optional maps)
# =============================================================================


# =============================================================================
# 1. PACKAGES, FUNCTIONS AND PARAMETERS ----
# =============================================================================

pacman::p_load('tidyverse', 'lubridate', 'readr', 'sf', 'adehabitatHR',
               'track2KBA', 'ggplot2', 'rnaturalearth', 'rnaturalearthdata')

source('functions/phenology_dates_functions.R')
source('parameters.R')


# =============================================================================
# 2. SETTINGS OF THIS SCRIPT ----
# =============================================================================

# Data to analyse (the same as in script 01):
#   'GLS'  example data: Cory's shearwater GLS tracks (input/)
#   'GPS'  example data: Audouin's gull GPS track (input/)
#   'full' full dataset of the paper (data/, not shared; created with
#          analyses/prepare_full_dataset.R). Results are saved in output/full/
dataset <- 'GLS'
if (dataset == 'full') {
  metadata_file  <- 'data/All_metadata.csv'
  out_dir        <- 'output/full'
} else {
  metadata_file  <- paste0('input/', dataset, '_metadata_example.csv')
  out_dir        <- 'output'
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Use the corrected classification (output of script 02)? If FALSE, the
# automated classification of script 01 is used.
use_corrected <- FALSE
tracks_file <- file.path(out_dir, ifelse(use_corrected,
                                         'Phenology_tracks_corrected.csv',
                                         'Phenology_tracks.csv'))

# OPTIONAL: name each wintering area with a map of areas (e.g. marine
# ecoregions). Folder or file of the polygons (NULL: not used) and column with
# the name of each area. The name is the area with most positions of the
# wintering stage (NA if no position falls inside any area).
areas_file <- NULL            # e.g. 'input/wintering_areas' (Marine Ecoregions of
                              #   the World, adapted by Virginia et al. 2025)
area_col   <- 'ECOREGION'

# Save a map with the KDEs and dates of each track?
save_plots <- TRUE
plots_dir  <- file.path(out_dir, 'plots', 'phenology')


# =============================================================================
# 3. READ DATA ----
# =============================================================================

metadata  <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(geo_trip = as.character(geo_trip))
tracks.ph <- read_csv(tracks_file, show_col_types = FALSE) %>%
  mutate(ID = as.character(ID))
ids <- unique(tracks.ph$ID)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
areas <- if (is.null(areas_file)) NULL else sf::read_sf(areas_file)
sf_use_s2(FALSE)
if (save_plots) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 4. LOOP OVER TRACKS ----
# =============================================================================

phen_list <- list()
failed    <- data.frame(ID = character(), error = character())

for (i in seq_along(ids)) {

  info     <- metadata %>% filter(geo_trip == ids[i])
  track.ph <- tracks.ph %>% filter(ID == ids[i])
  cat('\n[', i, '/', length(ids), '] Track', ids[i], '\n')

  if (nrow(info) != 1) {
    message('  Track ', ids[i], ' not found (or repeated) in the metadata')
    failed <- rbind(failed, data.frame(ID = ids[i], error = 'metadata'))
    next
  }

  res <- tryCatch(
    suppressMessages(
      phenology_track(track.ph = track.ph, info = info,
                      min_stage_length = min_stage_length,
                      kde_level = kde_level, adv = adv,
                      use_equinox = use_equinox,
                      world = world, make_plot = save_plots,
                      areas = areas, area_col = area_col)),
    error = function(e) {
      message('  ERROR in track ', ids[i], ': ', conditionMessage(e))
      failed <<- rbind(failed, data.frame(ID = ids[i],
                                          error = conditionMessage(e)))
      NULL
    })
  if (is.null(res)) next

  phen_list[[ids[i]]] <- res$dates
  if (save_plots && !is.null(res$map)) {
    ggsave(file.path(plots_dir, paste0(ids[i], '_phenology.png')), res$map,
           width = 15, height = 18, units = 'cm', dpi = 150)
  }
}


# =============================================================================
# 5. SAVE RESULTS ----
# =============================================================================

phen_all <- bind_rows(phen_list)
write_csv(phen_all, file.path(out_dir, 'Phenology.csv'))
if (nrow(failed) > 0) write_csv(failed, file.path(out_dir, 'Failed_tracks_03.csv'))

cat('\nFinished:', nrow(phen_all), 'tracks with dates,', nrow(failed),
    'failed.\n')
