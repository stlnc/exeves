# Monthly P-E water cycle analysis
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/07b_implications_monthly.R
# Shows monthly P-E trajectories for all days, ExEvE days, and non-ExEvE days

library(data.table)
library(ggplot2)
library(ggpubr)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

axis_decimal <- function(x) sprintf("%.1f", x)

cat("Loading data...\n")
exeves    <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec      <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

# Prepare merge
prec_sub <- copy(prec[, .(grid_id, date, value)])
setnames(prec_sub, "value", "prec")
exeves_prec <- merge(exeves, prec_sub, by = c('grid_id', 'date'), all.x = TRUE)
setnames(exeves_prec, "value", "evap")
rm(prec_sub); gc()

#===============================================================================
# Helper: build monthly P-E data and classify changes
#===============================================================================
build_monthly_pe <- function(dt) {
  pe_sums <- unique(dt[, .(evap = sum(evap) / (SUB_PERIOD_YEARS * GRID_CELL_N),
                            prec = sum(prec) / (SUB_PERIOD_YEARS * GRID_CELL_N),
                            diff_pe = (sum(prec) - sum(evap)) / (SUB_PERIOD_YEARS * GRID_CELL_N)),
                       by = .(period, month(date))])
  pe_sums <- pe_sums[, .(prec, evap, diff_pe,
                          diff_prec = diff(prec),
                          diff_evap = diff(evap),
                          period), .(month)]
  pe_sums[, diff_diff_pe := diff_prec - diff_evap]
  pe_sums[, sum_diff_pe  := diff_prec + diff_evap]
  pe_sums[, mean_flux    := (prec + evap) / 2]
  
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
  
  # Ensure all condition levels exist
  needed_levels <- c('Wetter - Accelerated', 'Wetter - Decelerated',
                     'Drier - Accelerated', 'Drier - Decelerated')
  missing <- setdiff(needed_levels, levels(to_plot$Conditions))
  if (length(missing) > 0) {
    levels(to_plot$Conditions) <- c(levels(to_plot$Conditions), missing)
  }
  
  # Handle any months where conditions remain "Unknown" (e.g. edge case)
  if (any(to_plot$Conditions == "Unknown")) {
    to_plot[Conditions == "Unknown" & Period == "Up to 2001",
            Conditions := factor('Wetter - Decelerated')]
  }
  
  to_plot <- to_plot[order(month)]
  return(to_plot)
}

#===============================================================================
# Panel A: All days
#===============================================================================
cat("Panel A: All days (monthly)...\n")
to_plot <- build_monthly_pe(exeves_prec)

month_names <- data.frame(
  x = to_plot[Period == "Up to 2001", diff_pe],
  y = to_plot[Period == "Up to 2001", mean_flux],
  text = month.abb
)

gg_all <- ggplot(to_plot) +
  geom_point(aes(y = mean_flux, x = diff_pe, fill = Period, shape = Period),
             colour = "transparent", size = 2) +
  geom_line(aes(y = mean_flux, x = diff_pe, group = factor(month), col = Conditions),
            alpha = 0.5) +
  geom_text(data = month_names, aes(x, y, label = text),
            cex = 3.5, nudge_x = 2, nudge_y = 2, col = 'grey40') +
  scale_fill_manual(values = c('grey60', 'grey20')) +
  scale_color_manual(values = WATER_CYCLE_CHANGE_PALETTE) +
  scale_shape_manual(values = c(22, 21)) +
  xlab(expression(atop(P - E ~ "[mm/month]"))) +
  ylab(expression(atop((P + E) / 2 ~ " [mm/month]"))) +
  scale_y_continuous(labels = axis_decimal) +
  theme_linedraw() +
  theme(axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

#===============================================================================
# Panel B: ExEvE days
#===============================================================================
cat("Panel B: ExEvE days (monthly)...\n")
to_plot <- build_monthly_pe(exeves_prec[!is.na(event_80_95_id)])

month_names <- data.frame(
  x = to_plot[Period == "Up to 2001", diff_pe],
  y = to_plot[Period == "Up to 2001", mean_flux],
  text = month.abb
)

gg_event <- ggplot(to_plot) +
  geom_point(aes(y = mean_flux, x = diff_pe, fill = Period, shape = Period),
             colour = "transparent", size = 2) +
  geom_line(aes(y = mean_flux, x = diff_pe, group = month, col = Conditions),
            alpha = 0.5) +
  geom_text(data = month_names, aes(x, y, label = text),
            cex = 3.5, nudge_x = -0.5, nudge_y = 0.5, col = 'grey40') +
  scale_fill_manual(values = c('grey60', 'grey20')) +
  scale_color_manual(values = WATER_CYCLE_CHANGE_PALETTE[c(1, 3, 4, 2)]) +
  scale_shape_manual(values = c(22, 21)) +
  xlab(expression(atop(P - E ~ "[mm/month]"))) +
  ylab(expression(atop((P + E) / 2 ~ " [mm/month]"))) +
  scale_y_continuous(labels = axis_decimal) +
  theme_linedraw() +
  theme(axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

#===============================================================================
# Panel C: Non-ExEvE days
#===============================================================================
cat("Panel C: Non-ExEvE days (monthly)...\n")
to_plot <- build_monthly_pe(exeves_prec[is.na(event_80_95_id)])

month_names <- data.frame(
  x = to_plot[Period == "Up to 2001", diff_pe],
  y = to_plot[Period == "Up to 2001", mean_flux],
  text = month.abb
)

gg_not_event <- ggplot(to_plot) +
  geom_point(aes(y = mean_flux, x = diff_pe, fill = Period, shape = Period),
             colour = "transparent", size = 2) +
  geom_line(aes(y = mean_flux, x = diff_pe, group = month, col = Conditions),
            alpha = 0.5) +
  geom_text(data = month_names, aes(x, y, label = text),
            cex = 3.5, nudge_x = 1.5, nudge_y = 1.5, col = 'grey40') +
  scale_fill_manual(values = c('grey60', 'grey20')) +
  scale_color_manual(values = WATER_CYCLE_CHANGE_PALETTE[c(1, 3, 4, 2)]) +
  scale_shape_manual(values = c(22, 21)) +
  xlab(expression(atop(P - E ~ "[mm/month]"))) +
  ylab(expression(atop((P + E) / 2 ~ " [mm/month]"))) +
  scale_y_continuous(labels = axis_decimal) +
  theme_linedraw() +
  theme(axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

#===============================================================================
# Combine & save
#===============================================================================
cat("Saving figure...\n")
ggarrange(gg_all, gg_event, gg_not_event,
          ncol = 1, nrow = 3,
          labels = c("A", "B", "C"),
          legend = 'right', common.legend = TRUE)
ggsave(paste0(PATH_OUTPUT_FIGURES, "implications_monthly.png"), width = 8, height = 12)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "implications_monthly.png"), "\n")
rm(exeves, prec, evap_grid, exeves_prec); gc()
