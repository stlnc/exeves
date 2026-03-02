# P-E water cycle budget analysis with waffle plots
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/07a_implications.R
# Computes P-E for all days, ExEvE days, and non-ExEvE days across two periods

library(data.table)
library(ggplot2)
library(ggpubr)
library(waffle)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

axis_decimal <- function(x) sprintf("%.1f", x)

cat("Loading data...\n")
exeves    <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec      <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

# Keyed update-join: add prec column in-place (no full copy)
setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.value, on = .(grid_id, date)]
setnames(exeves, "value", "evap")
exeves_prec <- exeves
rm(prec); gc()

#===============================================================================
# Helper: compute P-E summaries and classify water cycle changes
#===============================================================================
build_pe_plot <- function(dt, evap_grid_dt) {
  pe_sums <- unique(dt[, .(evap = sum(evap) / (SUB_PERIOD_YEARS * DAYS_IN_YEAR),
                            prec = sum(prec) / (SUB_PERIOD_YEARS * DAYS_IN_YEAR),
                            diff_pe = (sum(prec) - sum(evap)) / (SUB_PERIOD_YEARS * DAYS_IN_YEAR)),
                       by = .(grid_id, period)])
  pe_sums <- pe_sums[, .(prec, evap, diff_pe,
                          diff_prec = diff(prec),
                          diff_evap = diff(evap),
                          period), .(grid_id)]
  pe_sums[, diff_diff_pe := diff_prec - diff_evap]
  pe_sums[, sum_diff_pe  := diff_prec + diff_evap]
  pe_sums[, mean_flux    := (prec + evap) / 2]
  pe_sums <- merge(evap_grid_dt, pe_sums, by = "grid_id")
  
  to_plot <- copy(pe_sums)
  to_plot[, Conditions := factor("Unknown")]
  levels(to_plot$Conditions) <- c('Wetter - Accelerated', 'Wetter - Decelerated',
                                   'Drier - Accelerated', 'Drier - Decelerated')
  to_plot[sum_diff_pe > 0 & diff_diff_pe > 0, Conditions := factor('Wetter - Accelerated')]
  to_plot[sum_diff_pe < 0 & diff_diff_pe > 0, Conditions := factor('Wetter - Decelerated')]
  to_plot[sum_diff_pe > 0 & diff_diff_pe < 0, Conditions := factor('Drier - Accelerated')]
  to_plot[sum_diff_pe < 0 & diff_diff_pe < 0, Conditions := factor('Drier - Decelerated')]
  
  levels(to_plot$period) <- c("Up to 2001", "After 2001")
  setnames(to_plot, "period", "Period")
  
  # Ensure all levels exist for plotting (fix if a category is absent)
  needed_levels <- c('Wetter - Accelerated', 'Wetter - Decelerated',
                     'Drier - Accelerated', 'Drier - Decelerated')
  missing <- setdiff(needed_levels, levels(to_plot$Conditions))
  if (length(missing) > 0) {
    levels(to_plot$Conditions) <- c(levels(to_plot$Conditions), missing)
  }
  
  return(to_plot)
}

build_waffle_grob <- function(to_plot) {
  to_plot_nested <- to_plot[, .N, Conditions]
  to_plot_nested$N <- round(100 * to_plot_nested$N / sum(to_plot_nested$N), 0)
  gg_waffle <- ggplot(to_plot_nested, aes(fill = Conditions, values = N)) +
    geom_waffle(alpha = 0.5, size = 1, colour = "white", n_rows = 10) +
    scale_fill_manual(values = WATER_CYCLE_CHANGE_PALETTE[c(1, 3, 2, 4)]) +
    xlab("") + ylab("") +
    theme_linedraw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  ggplotGrob(gg_waffle)
}

#===============================================================================
# Panel A: All days
#===============================================================================
cat("Panel A: All days...\n")
to_plot_all <- build_pe_plot(exeves_prec, evap_grid)

gg_all <- ggplot(to_plot_all) +
  geom_point(aes(y = mean_flux, x = diff_pe, fill = Period, shape = Period),
             colour = "transparent", size = 2) +
  geom_line(aes(y = mean_flux, x = diff_pe, group = grid_id, col = Conditions), alpha = 0.5) +
  scale_fill_manual(values = c('grey60', 'grey20')) +
  scale_color_manual(values = WATER_CYCLE_CHANGE_PALETTE) +
  scale_shape_manual(values = c(22, 21)) +
  xlab(expression(atop(P - E ~ "[mm/day]"))) +
  ylab(expression(atop((P + E) / 2 ~ " [mm/day]"))) +
  scale_y_continuous(labels = axis_decimal) +
  theme_linedraw() +
  theme(axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

sub_grob <- build_waffle_grob(to_plot_all)
xr <- diff(range(to_plot_all$diff_pe, na.rm = TRUE))
yr <- diff(range(to_plot_all$mean_flux, na.rm = TRUE))
gg_all <- gg_all +
  annotation_custom(sub_grob,
                    xmin = max(to_plot_all$diff_pe, na.rm = TRUE) - 0.35 * xr,
                    xmax = max(to_plot_all$diff_pe, na.rm = TRUE),
                    ymin = max(to_plot_all$mean_flux, na.rm = TRUE) - 0.35 * yr,
                    ymax = max(to_plot_all$mean_flux, na.rm = TRUE))

#===============================================================================
# Panel B: ExEvE days
#===============================================================================
cat("Panel B: ExEvE days...\n")
to_plot_event <- build_pe_plot(exeves_prec[!is.na(event_80_95_id)], evap_grid)

gg_event <- ggplot(to_plot_event) +
  geom_point(aes(y = mean_flux, x = diff_pe, fill = Period, shape = Period),
             colour = "transparent", size = 2) +
  geom_line(aes(y = mean_flux, x = diff_pe, group = grid_id, col = Conditions), alpha = 0.5) +
  scale_fill_manual(values = c('grey60', 'grey20')) +
  scale_color_manual(values = WATER_CYCLE_CHANGE_PALETTE[c(1, 3, 4, 2)]) +
  scale_shape_manual(values = c(22, 21)) +
  xlab(expression(atop(P - E ~ "[mm/day]"))) +
  ylab(expression(atop((P + E) / 2 ~ " [mm/day]"))) +
  scale_y_continuous(labels = axis_decimal) +
  theme_linedraw() +
  theme(axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

sub_grob <- build_waffle_grob(to_plot_event)
xr <- diff(range(to_plot_event$diff_pe, na.rm = TRUE))
yr <- diff(range(to_plot_event$mean_flux, na.rm = TRUE))
gg_event <- gg_event +
  annotation_custom(sub_grob,
                    xmin = max(to_plot_event$diff_pe, na.rm = TRUE) - 0.35 * xr,
                    xmax = max(to_plot_event$diff_pe, na.rm = TRUE),
                    ymin = min(to_plot_event$mean_flux, na.rm = TRUE),
                    ymax = min(to_plot_event$mean_flux, na.rm = TRUE) + 0.35 * yr)

#===============================================================================
# Panel C: Non-ExEvE days
#===============================================================================
cat("Panel C: Non-ExEvE days...\n")
to_plot_non <- build_pe_plot(exeves_prec[is.na(event_80_95_id)], evap_grid)

gg_not_event <- ggplot(to_plot_non) +
  geom_point(aes(y = mean_flux, x = diff_pe, fill = Period, shape = Period),
             colour = "transparent", size = 2) +
  geom_line(aes(y = mean_flux, x = diff_pe, group = grid_id, col = Conditions), alpha = 0.5) +
  scale_fill_manual(values = c('grey60', 'grey20')) +
  scale_color_manual(values = WATER_CYCLE_CHANGE_PALETTE[c(2, 4, 1, 3)]) +
  scale_shape_manual(values = c(22, 21)) +
  xlab(expression(atop(P - E ~ "[mm/day]"))) +
  ylab(expression(atop((P + E) / 2 ~ " [mm/day]"))) +
  scale_y_continuous(labels = axis_decimal) +
  theme_linedraw() +
  theme(axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

sub_grob <- build_waffle_grob(to_plot_non)
xr <- diff(range(to_plot_non$diff_pe, na.rm = TRUE))
yr <- diff(range(to_plot_non$mean_flux, na.rm = TRUE))
gg_not_event <- gg_not_event +
  annotation_custom(sub_grob,
                    xmin = max(to_plot_non$diff_pe, na.rm = TRUE) - 0.35 * xr,
                    xmax = max(to_plot_non$diff_pe, na.rm = TRUE),
                    ymin = max(to_plot_non$mean_flux, na.rm = TRUE) - 0.35 * yr,
                    ymax = max(to_plot_non$mean_flux, na.rm = TRUE))

#===============================================================================
# Combine & save
#===============================================================================
cat("Saving figure...\n")
ggarrange(gg_all, gg_event, gg_not_event,
          ncol = 1, nrow = 3,
          labels = c("A", "B", "C"),
          legend = 'right', common.legend = TRUE)
ggsave(paste0(PATH_OUTPUT_FIGURES, "implications.png"), width = 8, height = 12)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "implications.png"), "\n")
rm(exeves, prec, evap_grid, exeves_prec); gc()
