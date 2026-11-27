# =============================================================================
# CANDIDATE TRACKS FOR FIGURE 3 (contact sheets)
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Puts the phenology maps of the full dataset (made by
#              03_Phenology_dates.R with dataset = 'full') in contact sheets,
#              one or more pages per species, to choose the examples of
#              Figure 3. Each map is labelled with the track ID, the number of
#              wintering areas, the review status and the four main dates.
#              Tracks with all four dates and without flags are shown first.
#
# Input:  output/full/Phenology.csv
#         output/full/plots/phenology/<ID>_phenology.png
# Output: <in_dir>/plots/figure_candidates.pdf (not shared in GitHub)
#
# Run from the project root, after 03_Phenology_dates.R with dataset = 'full'.
# =============================================================================

pacman::p_load('tidyverse', 'png', 'grid')

in_dir   <- file.path('output', 'full')
maps_dir <- file.path(in_dir, 'plots', 'phenology')
out_file <- file.path(in_dir, 'plots', 'figure_candidates.pdf')
per_page <- 12     # maps per page (4 columns x 3 rows)

species_names <- c(CALBOR = "Cory's shearwater", CALDIO = "Scopoli's shearwater",
                   CALEDW = 'Cape Verde shearwater', PUFPUF = 'Manx shearwater',
                   PUFYEL = 'Yelkouan shearwater', HYDCAS = 'Band-rumped storm petrel',
                   HYDJAB = 'Cape Verde storm petrel')

phen <- read_csv(file.path(in_dir, 'Phenology.csv'), show_col_types = FALSE,
                 col_types = cols(ID = 'c')) %>%
  mutate(complete = !is.na(B1_dep_KDE) & !is.na(W_total_arr_KDE) &
           !is.na(W_total_dep_KDE) & !is.na(B2_arr_KDE),
         file = file.path(maps_dir, paste0(ID, '_phenology.png'))) %>%
  filter(file.exists(file)) %>%
  # tracks with all dates and without flags first
  arrange(Sp, desc(complete), Check != 'Ok', n_WA, ID)

fmt <- function(d) ifelse(is.na(d), '-', format(as.Date(d), '%d %b'))

pdf(out_file, width = 16, height = 13)
for (sp in unique(phen$Sp)) {
  d <- phen %>% filter(Sp == sp)
  pages <- split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / per_page))
  for (k in seq_along(pages)) {
    grid.newpage()
    title <- paste0(sp, ' - ', species_names[sp], ' (', nrow(d), ' tracks)',
                    if (length(pages) > 1) paste0(', page ', k, '/', length(pages)))
    pushViewport(viewport(layout = grid.layout(4, 4, heights = unit(c(0.6, 1, 1, 1),
                                                                   c('in', 'null', 'null', 'null')))))
    grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4),
              gp = gpar(fontsize = 16, fontface = 'bold'))
    for (j in seq_along(pages[[k]])) {
      r <- d[pages[[k]][j], ]
      row <- (j - 1) %/% 4 + 2; col <- (j - 1) %% 4 + 1
      pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
      grid.raster(readPNG(r$file), y = unit(0.45, 'npc'), height = unit(0.85, 'npc'))
      lab1 <- paste0(r$ID, ' | ', r$Colony, ' ', r$Year1, ' | WA: ', r$n_WA,
                     ' | ', r$Check)
      lab2 <- paste0('B dep ', fmt(r$B1_dep_KDE), ' - W arr ', fmt(r$W_total_arr_KDE),
                     ' - W dep ', fmt(r$W_total_dep_KDE), ' - B arr ', fmt(r$B2_arr_KDE))
      grid.text(lab1, y = unit(0.97, 'npc'),
                gp = gpar(fontsize = 9, fontface = 'bold',
                          col = ifelse(r$Check == 'Ok' & r$complete, 'black', 'firebrick')))
      grid.text(lab2, y = unit(0.92, 'npc'), gp = gpar(fontsize = 8))
      popViewport()
    }
    popViewport()
  }
}
invisible(dev.off())
cat('Contact sheets saved in', out_file, '\n')
