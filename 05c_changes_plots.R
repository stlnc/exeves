# Change visualisation: annual time series, monthly polar plot, spatial map
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/05c_changes_plots.R

library(data.table)
library(ggplot2)
library(ggpubr)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

# Helper
equal_breaks <- function(n = 3, s = 0.05, ...) {
  function(x) {
    d <- s * diff(range(x)) / (1 + 2 * s)
    seq(min(x) + d, max(x) - d, length = n)
  }
}

#===============================================================================
# DATA
#===============================================================================
evap <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
monthly_changes <- readRDS(file = paste0(PATH_OUTPUT, 'monthly_changes.rds'))
spatial_changes <- readRDS(file = paste0(PATH_OUTPUT, 'spatial_changes.rds'))
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

#===============================================================================
# ANNUAL TIME SERIES
#===============================================================================
cat("Creating annual time-series plots...\n")

to_plot <- evap[, .(`All days` = mean(value)), year(date)]
to_plot <- merge(to_plot, evap[!is.na(event_80_95_id), .(ExEvEs = mean(value)), year(date)], by = 'year')
to_plot <- merge(to_plot, evap[is.na(event_80_95_id), .(`Non-ExEvEs` = mean(value)), year(date)], by = 'year')
to_plot <- melt(to_plot, id.vars = 'year', variable.name = "Conditions")

gg_intensity <- ggplot(to_plot) +
  geom_line(aes(x = year, y = value, col = Conditions)) +
  geom_point(aes(x = year, y = value, col = Conditions)) +
  geom_vline(xintercept = year(END_PERIOD_1) + 1, col = 'grey60', linetype = 2) +
  xlab("Year") + ylab("Evaporation (mm/day)") +
  scale_color_manual(values = SUBDUED_PROF_PALETTE[c(2, 4, 1)]) +
  theme_linedraw()

to_plot <- evap[, .(`All days` = sum(value) / GRID_CELL_N), year(date)]
to_plot <- merge(to_plot, evap[!is.na(event_80_95_id), .(ExEvEs = sum(value) / GRID_CELL_N), year(date)], by = 'year')
to_plot <- merge(to_plot, evap[is.na(event_80_95_id), .(`Non-ExEvEs` = sum(value) / GRID_CELL_N), year(date)], by = 'year')
to_plot <- melt(to_plot, id.vars = 'year', variable.name = "Conditions")

gg_severity <- ggplot(to_plot) +
  geom_line(aes(x = year, y = value, col = Conditions)) +
  geom_point(aes(x = year, y = value, col = Conditions)) +
  geom_vline(xintercept = year(END_PERIOD_1) + 1, col = 'grey60', linetype = 2) +
  xlab("Year") + ylab("Evaporation (mm/year)") +
  scale_color_manual(values = SUBDUED_PROF_PALETTE[c(2, 4, 1)]) +
  theme_linedraw()

#===============================================================================
# MONTHLY POLAR PLOT (evaporation)
#===============================================================================
cat("Creating monthly ratio plots...\n")

dummy <- melt(monthly_changes, id.vars = c('grid_id', 'month', 'period', 'conditions'))
dummy <- dcast(dummy, grid_id + month + conditions + variable ~ period)
dummy[, total_value_up_to_2001 := sum(up_to_2001, na.rm = TRUE), .(grid_id, variable, month)]
dummy[, total_value_after_2001 := sum(after_2001, na.rm = TRUE), .(grid_id, variable, month)]
dummy[, ratio_exeves := after_2001 / up_to_2001]
dummy[, ratio_total  := total_value_after_2001 / total_value_up_to_2001]
dummy <- dummy[conditions == "ExEvE" & variable %in% c("evap"),
               .(grid_id, month, variable, ratio_total, ratio_exeves)]

to_plot_1 <- dummy[, .(median = median(ratio_total, na.rm = TRUE),
                        q95 = quantile(ratio_total, 0.99, na.rm = TRUE),
                        q05 = quantile(ratio_total, 0.01, na.rm = TRUE),
                        Conditions = factor('All days')), .(month, variable)]
to_plot_2 <- dummy[, .(median = median(ratio_exeves, na.rm = TRUE),
                        q95 = quantile(ratio_exeves, 0.99, na.rm = TRUE),
                        q05 = quantile(ratio_exeves, 0.01, na.rm = TRUE),
                        Conditions = factor('ExEvEs')), .(month, variable)]
to_plot <- rbind(to_plot_1, to_plot_2)
to_plot <- melt(to_plot, id.vars = c('month', 'variable', 'Conditions'), variable.name = 'stat')

monthly_plot_sum <- ggplot() +
  geom_hline(aes(yintercept = y), data.frame(y = c(1:7)), color = "lightgrey") +
  geom_hline(aes(yintercept = y), data.frame(y = c(0, 8)), color = '#536878', alpha = .9) +
  geom_col(data = to_plot[stat == 'median' & Conditions == 'ExEvEs'],
           aes(x = month, y = value, fill = value),
           position = "dodge2", show.legend = TRUE, alpha = .8) +
  geom_point(data = to_plot[stat == 'median' & Conditions == 'All days'],
             aes(x = month, y = value), size = 1.5, color = '#536878') +
  geom_line(data = to_plot[Conditions == 'All days'],
            aes(x = month, y = value), color = '#536878') +
  geom_segment(data = to_plot[stat == 'median' & Conditions == 'ExEvEs'],
               aes(x = month, y = 0, xend = month, yend = 8),
               color = '#536878', linetype = 'dotted') +
  coord_polar() +
  scale_y_continuous(limits = c(-1.5, 8.5)) +
  scale_fill_gradientn("Ratio", colours = c('grey97', colset_subdued_prof[2])) +
  guides(fill = guide_colorsteps(barwidth = 9, barheight = .5,
                                  title.position = "top", title.hjust = .5)) +
  theme(axis.title = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = colset_subdued_prof[3], size = 10, vjust = 0.1),
        legend.position = "bottom",
        legend.margin = margin(-15, 0, 10, 0),
        legend.box.margin = margin(-5, 0, 5, 0),
        text = element_text(color = colset_subdued_prof[3]),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.grid = element_blank(),
        plot.margin = margin(b = 0, l = -2, unit = "cm"))

#===============================================================================
# SPATIAL PLOT
#===============================================================================
cat("Creating spatial plot...\n")
setnames(spatial_changes, "ratio", "Ratio")

lon_range <- evap_grid[, range(lon)]
lat_range <- evap_grid[, range(lat)]

spatial_plot <- ggplot(spatial_changes[period == "up_to_2001" &
                                       variable == "Evaporation (ExEvEs)"]) +
  geom_tile(aes(lon, lat, fill = Ratio)) +
  scale_fill_gradientn("Ratio", colours = c('grey97', SUBDUED_PROF_PALETTE[4])) +
  xlab('') + ylab('') +
  theme_linedraw() +
  theme(text = element_text(color = SUBDUED_PROF_PALETTE[3]),
        legend.direction = "horizontal", legend.position = "bottom",
        legend.margin = margin(-15, 0, 10, 0),
        legend.box.margin = margin(-5, 0, 5, 0),
        plot.margin = margin(t = 1, b = 0, r = 1, unit = "cm")) +
  guides(fill = guide_colorsteps(title.position = "top", title.hjust = 0.5,
                                  barwidth = 9, barheight = .5))

#===============================================================================
# COMBINED FIGURE
#===============================================================================
gg_1 <- ggarrange(gg_severity, gg_intensity, nrow = 1, labels = c("A", "B"),
                   legend = 'right', common.legend = TRUE)
gg_2 <- ggarrange(spatial_plot, monthly_plot_sum, nrow = 1,
                   labels = c("C", "D"), widths = c(1, 1),
                   legend = 'bottom', common.legend = FALSE)

ggarrange(gg_1, gg_2, nrow = 2) + bgcolor("white")
ggsave(paste0(PATH_OUTPUT_FIGURES, "exeve_changes.pdf"), width = 10, height = 8, bg = "white")
ggsave(paste0(PATH_OUTPUT_FIGURES, "exeve_changes.png"), width = 10, height = 8, bg = "white", dpi = 300)

cat("Change plots saved.\n")
rm(evap, monthly_changes, spatial_changes); gc()
