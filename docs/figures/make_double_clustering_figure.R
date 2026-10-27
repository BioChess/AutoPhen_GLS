# =============================================================================
# Figure: double hierarchical clustering of segments (docs/Classification_rules.md)
# =============================================================================
# Run from the project root. Output: docs/figures/double_clustering.png

pacman::p_load('tidyverse', 'patchwork', 'gridExtra', 'grid', 'ragg')
Sys.setlocale('LC_TIME', 'C')   # month names in English
invisible(suppressWarnings(Sys.setlocale('LC_CTYPE', 'C.UTF-8')))  # symbols (Linux)
font <- 'Liberation Sans'       # use 'Arial' on Windows

col_B <- '#E69500'   # breeding (orange)
col_W <- '#4682B4'   # wintering (blue)

# --- 1. Distance to the colony of an example track ----------------------------
# Positions of the example track after step 1 (output of script 01): x_km and
# y_km are projected coordinates centred on the colony
id  <- '11001001'
trk <- read_csv('output/Phenology_tracks.csv', show_col_types = FALSE,
                col_types = cols(ID = 'c')) %>%
  filter(ID == id) %>%
  mutate(dist_km = sqrt(x_km^2 + y_km^2)) %>%
  arrange(date)

p_dist <- ggplot(trk, aes(date, dist_km, colour = dist_km)) +
  geom_line(linewidth = 0.9) +
  scale_colour_gradient(low = col_B, high = col_W, guide = 'none') +
  scale_x_datetime(date_breaks = '1 month', date_labels = '%b') +
  labs(x = 'Date', y = 'Distance to the colony (km)') +
  theme_bw(base_size = 11, base_family = font) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 2. Thresholds of the movement criteria -----------------------------------
tab <- data.frame(
  ` ` = c('Stationary', 'Migration'),
  Straightness = c('≤ 0.5', '> 0.5'),
  `Mean vector of\nturning angles` = c('≤ 0.5', '> 0.5'),
  Sinuosity = c('≥ 1', '< 1'),
  `Speed (mean\nstep length)` = c('≤ 1.5 × median', '> 1.5 × median'),
  check.names = FALSE)
tt <- ttheme_minimal(
  core    = list(fg_params = list(fontsize = 11, fontfamily = font),
                 bg_params = list(fill = 'grey95', col = 'black')),
  colhead = list(fg_params = list(fontsize = 11, fontface = 'bold', fontfamily = font),
                 bg_params = list(fill = 'grey75', col = 'black')),
  rowhead = list(fg_params = list(fontsize = 11, fontface = 'bold', fontfamily = font)))
tg <- tableGrob(tab, rows = NULL, theme = tt)
# First column (row names) in bold with a grey background
for (i in which(tg$layout$l == 1 & tg$layout$t > 1)) {
  if (grepl('core-fg', tg$layout$name[i])) tg$grobs[[i]]$gp <- gpar(fontsize = 11, fontface = 'bold', fontfamily = font)
  if (grepl('core-bg', tg$layout$name[i])) tg$grobs[[i]]$gp <- gpar(fill = 'grey75', col = 'black')
}
p_tab <- wrap_elements(full = tg) + theme(plot.margin = margin(5, 5, 5, 15))

# --- 3. Diagram of the two clusterings ----------------------------------------
box <- function(x, y, w, h, label, col = 'black', size = 6, border = 'black',
                fill = 'white') {
  list(annotate('rect', xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2,
                ymax = y + h / 2, fill = fill, colour = border),
       annotate('text', x = x, y = y, label = label, colour = col,
                size = size, fontface = 'bold', lineheight = 0.9,
                family = font))
}
seg <- function(x, y, xend, yend) annotate('segment', x = x, y = y,
                                           xend = xend, yend = yend)
xB <- 2; xW <- 5; xM <- 8
p_diag <- ggplot() +
  # Movement mode: Stationary (B, W) vs Migration
  seg(6.5, 8.6, 6.5, 8.1) + seg(3.5, 8.1, 8, 8.1) +
  seg(3.5, 8.1, 3.5, 7.4) + seg(xB, 7.4, xW, 7.4) +
  seg(xB, 7.4, xB, 6.4) + seg(xW, 7.4, xW, 6.4) + seg(xM, 8.1, xM, 6.4) +
  # Distance: Breeding vs Non-breeding (W, mig)
  seg(xB, 5.2, xB, 3.6) + seg(xW, 5.2, xW, 4.4) + seg(xM, 5.2, xM, 4.4) +
  seg(xW, 4.4, xM, 4.4) + seg(6.5, 4.4, 6.5, 3.6) + seg(xB, 3.6, 6.5, 3.6) +
  seg(4.25, 3.6, 4.25, 3.0) +
  box(6.5, 9.4, 4.2, 1.6, 'Stationary vs\nmigration') +
  box(xB, 5.8, 2.6, 1.0, 'Breeding', col = col_B, border = 'darkgreen',
      fill = 'grey95', size = 5.5) +
  box(xW, 5.8, 2.6, 1.0, 'Wintering', col = col_W, border = 'darkgreen',
      fill = 'grey95', size = 5.5) +
  box(xM, 5.8, 2.6, 1.0, 'Migration', border = 'darkgreen', fill = 'grey95',
      size = 5.5) +
  box(4.25, 2.2, 4.2, 1.6, 'Breeding vs\nnon-breeding') +
  coord_cartesian(xlim = c(0.5, 9.5), ylim = c(1.2, 10.4), expand = FALSE) +
  theme_void()

# --- 4. Combine and save ------------------------------------------------------
fig <- (p_tab / p_dist + plot_layout(heights = c(1, 2.2))) | p_diag
fig <- fig + plot_layout(widths = c(1, 1.15))
ggsave('docs/figures/double_clustering.png', fig, width = 30, height = 16,
       units = 'cm', dpi = 200, bg = 'white', device = ragg::agg_png)
