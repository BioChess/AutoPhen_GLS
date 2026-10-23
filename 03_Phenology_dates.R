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

# Use the same example data as in script 01: 'GLS' or 'GPS'
example <- 'GLS'
metadata_file <- paste0('input/', example, '_metadata_example.csv')

# Use the corrected classification (output of script 02)? If FALSE, the
# automated classification of script 01 is used.
use_corrected <- FALSE
tracks_file <- ifelse(use_corrected, 'output/Phenology_tracks_corrected.csv',
                      'output/Phenology_tracks.csv')

# Save a map with the KDEs and dates of each track?
save_plots <- TRUE
plots_dir  <- 'output/plots/phenology'


# =============================================================================
# 3. READ DATA ----
# =============================================================================

metadata  <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(geo_trip = as.character(geo_trip))
tracks.ph <- read_csv(tracks_file, show_col_types = FALSE) %>%
  mutate(ID = as.character(ID))
ids <- unique(tracks.ph$ID)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
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
                      world = world, make_plot = save_plots)),
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
write_csv(phen_all, 'output/Phenology.csv')
if (nrow(failed) > 0) write_csv(failed, 'output/Failed_tracks_03.csv')

cat('\nFinished:', nrow(phen_all), 'tracks with dates,', nrow(failed),
    'failed.\n')
