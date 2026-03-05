# Period-split day-of-extreme analysis
# Extends 07c by splitting Panels A and B by period (up_to_2001 vs after_2001)
# to reveal whether the timing of extreme evap / prec peak within events has shifted
#
# Output: figures/wet_days_by_period.png

library(data.table)
library(ggplot2)
library(ggpubr)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec   <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))

# Add prec via keyed update-join
setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.value, on = .(grid_id, date)]
setnames(exeves, "value", "evap")
rm(prec); gc()

# Subset to ExEvE events
exeves_prec <- exeves[!is.na(event_80_95_id)]
exeves_prec[!is.na(extreme_id), extreme_evap := TRUE]

#===============================================================================
# EVENT-DAY METRICS (same as 07c)
#===============================================================================
cat("Computing event-day metrics...\n")

exeves_prec[, event_day := seq_len(.N), by = .(event_80_95_id, grid_id)]
exeves_prec[, extreme_day := seq_len(.N), by = .(extreme_id, grid_id)]
exeves_prec[, event_duration := .N, by = .(event_80_95_id, grid_id)]
exeves_prec[is.na(extreme_id), extreme_day := NA]
exeves_prec[extreme_day == 1, event_day_of_extreme := event_day]

# Day of precipitation maximum within each event
prec_max_idx <- exeves_prec[, .I[which.max(prec)], .(event_80_95_id, grid_id)]$V1
prec_max <- exeves_prec[prec_max_idx, .(date, grid_id, event_80_95_id, prec_max = prec)]
exeves_prec <- merge(exeves_prec, prec_max, by = c('grid_id', 'event_80_95_id', 'date'), all.x = TRUE)
exeves_prec[!is.na(prec_max), event_day_of_prec_max := event_day]

exeves_prec[, event_day_of_first_extreme := suppressWarnings(min(event_day_of_extreme, na.rm = TRUE)),
            .(grid_id, event_80_95_id)]
exeves_prec[, event_day_of_first_prec_max := suppressWarnings(min(event_day_of_prec_max, na.rm = TRUE)),
            .(grid_id, event_80_95_id)]
exeves_prec[is.infinite(event_day_of_first_extreme), event_day_of_first_extreme := NA_real_]
exeves_prec[is.infinite(event_day_of_first_prec_max), event_day_of_first_prec_max := NA_real_]

# Duration selection (same as 07c)
MIN_VIZ_DURATION <- 6L
dur_counts <- exeves_prec[, .(N = uniqueN(paste(grid_id, event_80_95_id))), by = event_duration]
dur_above <- dur_counts[event_duration >= MIN_VIZ_DURATION]
if (nrow(dur_above) > 0) {
  MODAL_DURATION <- dur_above[which.max(N), event_duration]
} else {
  MODAL_DURATION <- dur_counts[which.max(N), event_duration]
}
cat("  Visualisation duration:", MODAL_DURATION, "days\n")

#===============================================================================
# Panel A: Day-of-extreme density — SPLIT BY PERIOD
#===============================================================================
cat("Panel A: Day-of-extreme density by period...\n")

to_plot <- unique(exeves_prec[, .(grid_id, event_80_95_id, event_duration,
                                   event_day_of_first_extreme,
                                   event_day_of_first_prec_max,
                                   month = month(date),
                                   period)])

# Relabel periods for display
to_plot[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
to_plot[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]

gg_extreme_day <- ggplot(to_plot[event_duration == MODAL_DURATION]) +
  geom_vline(xintercept = 1, col = 'black') +
  geom_hline(yintercept = 0, col = 'black') +
  geom_density(aes(x = event_day_of_first_prec_max, group = interaction(event_duration, Period),
                   col = Period), stat = "count", linetype = 2) +
  geom_density(aes(x = event_day_of_first_extreme, group = interaction(event_duration, Period),
                   col = Period), stat = "count") +
  facet_wrap(~ month, scales = 'free_y', ncol = 4) +
  scale_color_manual(values = c("Up to 2001" = PALETTES$subdued_prof[2],
                                 "After 2001" = PALETTES$subdued_prof[4])) +
  xlab("Day of event") +
  ylab("Number of events") +
  labs(subtitle = "Solid = first extreme evap day; Dashed = precipitation peak day") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

#===============================================================================
# Panel B: Mean evap and prec by event day — SPLIT BY PERIOD
#===============================================================================
cat("Panel B: Mean fluxes by event day and period...\n")

to_plot_mean <- exeves_prec[, .(Evaporation = mean(evap, na.rm = TRUE),
                                 Precipitation = mean(prec, na.rm = TRUE)),
                             .(event_duration, event_day, month = month(date), period)]
to_plot_mean <- melt(to_plot_mean, id.vars = c("event_duration", "event_day", "month", "period"))
setnames(to_plot_mean, "variable", "Variable")
to_plot_mean[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
to_plot_mean[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]

gg_mean_day <- ggplot(to_plot_mean[event_duration == MODAL_DURATION]) +
  geom_vline(xintercept = 1, col = 'black') +
  geom_hline(yintercept = 0, col = 'black') +
  geom_line(aes(x = event_day, y = value, col = Variable, linetype = Period)) +
  facet_wrap(~ month, scales = 'free_y', ncol = 4) +
  xlab("Day of event") +
  ylab("Mean (mm/day)") +
  scale_color_manual(values = c("Evaporation" = PALETTES$subdued_prof[4],
                                 "Precipitation" = PALETTES$subdued_prof[2])) +
  scale_linetype_manual(values = c("Up to 2001" = 1, "After 2001" = 2)) +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

#===============================================================================
# Panel C: Lag between prec peak and extreme evap — by period
#===============================================================================
cat("Panel C: Lag distribution by period...\n")

lag_dt <- unique(exeves_prec[event_duration == MODAL_DURATION,
                              .(grid_id, event_80_95_id,
                                event_day_of_first_extreme,
                                event_day_of_first_prec_max,
                                month = month(date),
                                period)])
lag_dt[, lag := event_day_of_first_prec_max - event_day_of_first_extreme]
lag_dt <- lag_dt[!is.na(lag)]
lag_dt[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
lag_dt[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]

gg_lag <- ggplot(lag_dt) +
  geom_bar(aes(x = lag, fill = Period), position = "dodge") +
  facet_wrap(~ month, scales = 'free_y', ncol = 4) +
  scale_fill_manual(values = c("Up to 2001" = PALETTES$subdued_prof[2],
                                "After 2001" = PALETTES$subdued_prof[4])) +
  geom_vline(xintercept = 0, linetype = 2, col = 'grey40') +
  xlab("Lag (prec peak day − extreme evap day)") +
  ylab("Number of events") +
  labs(subtitle = "Negative = prec peaks before evap extreme; Positive = prec peaks after") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

#===============================================================================
# Combine & save
#===============================================================================
cat("Saving figure...\n")
ggarrange(gg_extreme_day, gg_mean_day, gg_lag,
          ncol = 1, labels = c("A", "B", "C"),
          legend = 'bottom')
ggsave(paste0(PATH_OUTPUT_FIGURES, "wet_days_by_period.png"), width = 10, height = 14, dpi = 300)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "wet_days_by_period.png"), "\n")
rm(exeves, exeves_prec); gc()
