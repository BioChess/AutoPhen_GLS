# =============================================================================
# SENSITIVITY TESTS 04: FIGURES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Two kinds of figures:
#   (1) One panel per parameter, with the values of each track: agreement of
#       the classification (fig_sens_class) and mean error of the dates
#       (fig_sens_dates). The value of the main analysis is in colour.
#   (2) Summary of all tests (fig_sens_summary): mean paired difference with
#       the main analysis (variant - default) and bootstrap 95% CI, for the
#       agreement of the classification and for the dates within 3 days.
#
# Input:  output/sensitivity/_results/ (from 03_evaluate_variants.R)
# Output: output/sensitivity/_results/fig_sens_*.png and .pdf
# =============================================================================

pacman::p_load('tidyverse', 'patchwork', 'scales')
source('analyses/sensitivity/sensitivity_variants.R')

res_dir   <- file.path('output', 'sensitivity', '_results')
per_track <- read_csv(file.path(res_dir, 'sens_per_track.csv'),
                      show_col_types = FALSE, col_types = cols(ID = 'c'))
summ      <- read_csv(file.path(res_dir, 'sens_summary.csv'), show_col_types = FALSE)

col_default <- '#2a6f97'   # main analysis
col_variant <- 'grey55'    # other values

theme_sens <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = 'grey95', colour = NA),
        strip.text = element_text(face = 'bold'),
        legend.position = 'none')

save_fig <- function(p, name, width, height) {
  ggsave(file.path(res_dir, paste0(name, '.png')), p, width = width,
         height = height, units = 'cm', dpi = 300)
  ggsave(file.path(res_dir, paste0(name, '.pdf')), p, width = width,
         height = height, units = 'cm')
}


# =============================================================================
# 1. ONE PANEL PER PARAMETER ----
# =============================================================================
# The results of 'default' are repeated in each panel, with the default value
# of that parameter, so each panel compares the values of one parameter.

params <- per_track %>% filter(variant != 'default') %>%
  distinct(parameter, default, step)
default_rows <- per_track %>% filter(variant == 'default') %>% select(-parameter, -default, -value, -step)

panel_data <- bind_rows(
  per_track %>% filter(variant != 'default') %>% mutate(is_default = FALSE),
  crossing(params, default_rows) %>% mutate(value = default, is_default = TRUE)
)

# Order of the values inside each panel: numeric when possible
value_order <- function(d) {
  num <- suppressWarnings(as.numeric(d$value))
  ord <- if (any(!is.na(num))) order(num, na.last = TRUE) else order(!d$is_default, d$value)
  unique(d$value[ord])
}
# x values are 'parameter|value' (the same value can be in several panels);
# only the value is shown in the axis
panel_data <- panel_data %>%
  group_by(parameter) %>%
  mutate(x = paste(parameter, value, sep = '|'),
         x_order = match(value, value_order(pick(value, is_default)))) %>%
  ungroup() %>%
  arrange(parameter, x_order) %>%
  mutate(x = factor(x, levels = unique(x)),
         parameter = factor(parameter, levels = unique(
           map_chr(variants[-1], 'parameter'))))

panel_plot <- function(d, y, ylab, ytrans = 'identity', ybreaks = waiver()) {
  ggplot(d, aes(x, .data[[y]], fill = is_default, colour = is_default)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.25, linewidth = 0.4) +
    geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.5) +
    scale_fill_manual(values = c(`TRUE` = col_default, `FALSE` = col_variant)) +
    scale_colour_manual(values = c(`TRUE` = col_default, `FALSE` = col_variant)) +
    scale_x_discrete(labels = function(x) str_wrap(sub('.*[|]', '', x), 9)) +
    scale_y_continuous(trans = ytrans, breaks = ybreaks) +
    facet_wrap(~ parameter, scales = 'free_x', ncol = 4,
               labeller = label_wrap_gen(25)) +
    labs(x = NULL, y = ylab) +
    theme_sens
}

# Classification: only the parameters that change the classification
p_class <- panel_plot(panel_data %>% filter(step == '01'),
                      'class_agreement', 'Agreement with the reference classification\n(proportion of positions per track)')
save_fig(p_class, 'fig_sens_class', 18, 14)

# Dates: all parameters
p_dates <- panel_plot(panel_data, 'mean_abs_error',
                      'Mean absolute difference with the reference dates\n(days per track)',
                      ytrans = pseudo_log_trans(sigma = 1, base = 10),
                      ybreaks = c(0, 1, 2, 3, 5, 10, 20, 50, 100))
save_fig(p_dates, 'fig_sens_dates', 18, 18)


# =============================================================================
# 2. SUMMARY OF ALL TESTS ----
# =============================================================================

sum_long <- summ %>%
  filter(variant != 'default') %>%
  transmute(variant, group, parameter, value,
            d_class, d_class_low, d_class_high,
            d_within3, d_within3_low, d_within3_high, step) %>%
  pivot_longer(c(d_class, d_within3), names_to = 'metric', values_to = 'diff') %>%
  mutate(low  = ifelse(metric == 'd_class', d_class_low, d_within3_low),
         high = ifelse(metric == 'd_class', d_class_high, d_within3_high),
         # the classification does not change in the step '03' variants
         diff = ifelse(metric == 'd_class' & step == '03', NA, diff),
         metric = recode(metric,
                         d_class   = 'Agreement of the classification',
                         d_within3 = 'Dates within 3 days of the reference'),
         parameter = factor(parameter, levels = unique(
           map_chr(variants[-1], 'parameter'))),
         label = value)

p_sum <- ggplot(sum_long, aes(diff, label, colour = group)) +
  geom_vline(xintercept = 0, colour = 'grey40', linewidth = 0.4) +
  geom_errorbar(aes(xmin = low, xmax = high), width = 0, linewidth = 0.6,
                orientation = 'y',
                 na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_colour_manual(values = c(Requested = col_default, Additional = 'grey35'),
                      name = NULL) +
  facet_grid(parameter ~ metric, scales = 'free', space = 'free_y',
             switch = 'y', labeller = label_wrap_gen(22)) +
  labs(x = 'Difference with the main analysis\n(percentage points; mean and bootstrap 95% CI)',
       y = NULL) +
  theme_sens +
  theme(strip.placement = 'outside',
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = 'plain'),
        panel.grid.major.y = element_blank(),
        panel.spacing.x = unit(0.6, 'cm'),
        legend.position = 'bottom')
save_fig(p_sum, 'fig_sens_summary', 17, 22)
