# =============================================================================
# Panel of the figure of the double hierarchical clustering: distance to the
# colony of an example track (docs/Classification_rules.md)
# =============================================================================
# Run from the project root. Output: docs/figures/colony_distance_example.png
# and .svg (the rest of the figure is assembled by hand)

pacman::p_load('tidyverse', 'ragg')
Sys.setlocale('LC_TIME', 'C')   # month names in English
invisible(suppressWarnings(Sys.setlocale('LC_CTYPE', 'C.UTF-8')))  # symbols (Linux)
font <- 'Liberation Sans'       # use 'Arial' on Windows

col_B <- '#E69500'   # breeding (orange)
col_W <- '#4682B4'   # wintering (blue)

# --- 1. Distance to the colony of an example track ----------------------------
# Positions of the example track after step 1 (output of script 01): x_km and
# y_km are projected coordinates centred on the colony. As the figure is only
# illustrative, the distance is smoothed with a running median (9 positions,
# 4.5 days) and a running mean (3 positions).
id  <- '19001001'
trk <- read_csv('output/Phenology_tracks.csv', show_col_types = FALSE,
                col_types = cols(ID = 'c')) %>%
  filter(ID == id) %>%
  arrange(date) %>%
  mutate(dist_km = sqrt(x_km^2 + y_km^2),
         dist_km = runmed(dist_km, k = 9, endrule = 'median'),
         dist_km = as.numeric(stats::filter(dist_km, rep(1 / 3, 3), sides = 2)),
         dist_km = coalesce(dist_km, sqrt(x_km^2 + y_km^2)))

p_dist <- ggplot(trk, aes(date, dist_km, colour = dist_km)) +
  geom_line(linewidth = 1, lineend = 'round', linejoin = 'round') +
  scale_colour_gradient(low = col_B, high = col_W, guide = 'none') +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
  labs(x = 'Date', y = 'Distance to the colony (km)') +
  theme_bw(base_size = 11, base_family = font) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 2. Save -------------------------------------------------------------------
ggsave('docs/figures/colony_distance_example.png', p_dist, width = 14, height = 10,
       units = 'cm', dpi = 300, bg = 'white', device = ragg::agg_png)
ggsave('docs/figures/colony_distance_example.svg', p_dist, width = 14, height = 10,
       units = 'cm', bg = 'white', device = grDevices::svg)
