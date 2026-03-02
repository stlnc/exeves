# Drivers plot: evaporation vs precipitation under ExEvE / non-ExEvE conditions
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/06b_drivers_plot.R
# Uses precipitation as the available co-variable (radiation/temp not available)

library(data.table)
library(ggplot2)
library(ggpubr)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

cat("Loading drivers data...\n")
exeves_drivers <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_exeves_drivers.rds'))

#===============================================================================
# SEASONAL SCATTER: evaporation vs precipitation
#===============================================================================
to_plot <- exeves_drivers[month(date) %in% c(3, 6, 9, 12),
                          .(evap = mean(evap), prec = mean(prec)),
                          .(grid_id, month(date), conditions)]
to_plot[, month := month(month, label = TRUE)]

gg_prec <- ggplot(to_plot) +
  geom_point(aes(x = evap, y = prec, col = conditions), alpha = 0.5) +
  facet_wrap(~month, scales = 'free', ncol = 4) +
  ylab("Precip. (mm/day)") +
  xlab("Evaporation (mm/day)") +
  scale_color_manual(values = colset_subdued_prof[c(4, 2)]) +
  guides(col = guide_legend(title = "Conditions")) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 12),
        strip.background = element_rect(fill = 'grey30'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(paste0(PATH_OUTPUT_FIGURES, "drivers.png"), gg_prec, width = 10, height = 4)
cat("Saved drivers.png\n")

#===============================================================================
# STANDARDISED SCATTER: std evaporation vs precipitation
#===============================================================================
to_plot_std <- exeves_drivers[, .(evap = mean(std_value), prec = mean(prec)),
                               .(grid_id, month(date), conditions)]

gg_evap_prec_std <- ggplot(to_plot_std) +
  geom_point(aes(x = evap, y = prec, col = conditions), alpha = 0.7) +
  geom_hline(yintercept = 0, col = colset_subdued_prof[3]) +
  geom_vline(xintercept = 0, col = colset_subdued_prof[3]) +
  facet_wrap(~month, scales = 'free') +
  xlab("Evaporation (z-score)") +
  ylab("Precipitation (mm/day)") +
  scale_color_manual(values = colset_subdued_prof[c(4, 2)]) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -0.5,
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjust = -0.5,
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.background = element_rect(fill = 'grey30'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(paste0(PATH_OUTPUT_FIGURES, "std_evap_prec.png"), gg_evap_prec_std,
       width = 8, height = 9)
cat("Saved std_evap_prec.png\n")

rm(exeves_drivers); gc()
cat("Drivers plot analysis complete.\n")
