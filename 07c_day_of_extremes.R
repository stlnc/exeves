# Day-of-extreme within events analysis
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/07c_day_of_extremes.R
# Panel A: Density of first extreme / first prec-max day within events
# Panel B: Mean evap and prec by event day
# Panel C: Wet/Dry ratio by period and month

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

# Add prec to full exeves via keyed update-join (keep for Panel C)
setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.value, on = .(grid_id, date)]
setnames(exeves, "value", "evap")
rm(prec); gc()

# For Panels A&B: subset to ExEvE events only
exeves_prec <- exeves[!is.na(event_80_95_id)]

# Flag extreme evaporation days
exeves_prec[!is.na(extreme_id), extreme_evap := TRUE]

#===============================================================================
# 1. COMPUTE EVENT-DAY METRICS
#===============================================================================
cat("Computing event-day metrics...\n")

exeves_prec[, event_day := seq_len(.N), by = .(event_80_95_id, grid_id)]
exeves_prec[, extreme_day := seq_len(.N), by = .(extreme_id, grid_id)]
exeves_prec[, event_duration := .N, by = .(event_80_95_id, grid_id)]
exeves_prec[, extremes_per_event := .N, by = .(extreme_id, event_80_95_id, grid_id)]
exeves_prec[is.na(extreme_id), extreme_day := NA]
exeves_prec[is.na(extreme_id), extremes_per_event := NA]
exeves_prec[extreme_day == 1, event_day_of_extreme := event_day]

# Cumulative sums
exeves_prec[, cumsum_evap := cumsum(fifelse(is.na(evap), 0, evap)), .(grid_id, event_80_95_id)]
exeves_prec[, cumsum_prec := cumsum(fifelse(is.na(prec), 0, prec)), .(grid_id, event_80_95_id)]
exeves_prec[, cumsum_diff := cumsum_prec - cumsum_evap, .(grid_id, event_80_95_id)]

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

# Select the most common event duration for visualisation
dur_counts <- exeves_prec[, .(N = uniqueN(paste(grid_id, event_80_95_id))), by = event_duration]
MODAL_DURATION <- dur_counts[which.max(N), event_duration]
cat("  Most common event duration:", MODAL_DURATION, "days\n")

#===============================================================================
# Panel A: Density of first extreme / first prec-max day
#===============================================================================
cat("Panel A: Day-of-extreme density...\n")

to_plot <- unique(exeves_prec[, .(grid_id, event_80_95_id, event_duration,
                                   event_day_of_first_extreme,
                                   event_day_of_first_prec_max,
                                   month = month(date),
                                   period)])

gg_extreme_day <- ggplot(to_plot[event_duration == MODAL_DURATION]) +
  geom_vline(xintercept = 1, col = 'black') +
  geom_hline(yintercept = 0, col = 'black') +
  geom_density(aes(x = event_day_of_first_prec_max, group = event_duration), stat = "count",
               col = PALETTES$subdued_prof[2], linetype = 2) +
  geom_density(aes(x = event_day_of_first_extreme, group = event_duration), stat = "count",
               col = PALETTES$subdued_prof[4]) +
  facet_wrap(~ month, scales = 'free_y', ncol = 4) +
  xlab("Day of event") +
  ylab("Number of days") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

#===============================================================================
# Panel B: Mean evap and prec by event day
#===============================================================================
cat("Panel B: Mean fluxes by event day...\n")

to_plot_mean <- exeves_prec[, .(Evaporation = mean(evap, na.rm = TRUE), Precipitation = mean(prec, na.rm = TRUE)),
                             .(event_duration, event_day, month = month(date))]
to_plot_mean <- melt(to_plot_mean, id.vars = c("event_duration", "event_day", "month"))
setnames(to_plot_mean, "variable", "Variable")

gg_mean_day <- ggplot(to_plot_mean[event_duration == MODAL_DURATION]) +
  geom_vline(xintercept = 1, col = 'black') +
  geom_hline(yintercept = 0, col = 'black') +
  geom_line(aes(x = event_day, y = value, col = Variable, linetype = Variable)) +
  facet_wrap(~ month, scales = 'free_y', ncol = 4) +
  xlab("Day of event") +
  ylab("Mean (mm/day)") +
  scale_color_manual(values = PALETTES$subdued_prof[c(4, 2)]) +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

#===============================================================================
# Panel C: Wet/Dry day ratio
#===============================================================================
cat("Panel C: Wet/Dry day ratio...\n")

# Use the full exeves table (already has prec column from above)
exeves_all <- exeves

exeves_all[, prec_day := factor(fifelse(prec >= 1, "wet", "dry"), levels = c("wet", "dry"))]

# All days
all_prec_days <- exeves_all[!is.na(prec), .N / (GRID_CELL_N * 0.5 * PERIOD_LENGTH),
                             .(period, prec_day, month(date))]
setnames(all_prec_days, "V1", "value")
all_prec_days$conditions <- "All days"

# ExEvE days
exeves_prec_days <- exeves_all[!is.na(event_80_95_id) & !is.na(prec),
                                .N / (GRID_CELL_N * 0.5 * PERIOD_LENGTH),
                                .(period, prec_day, month(date))]
setnames(exeves_prec_days, "V1", "value")
exeves_prec_days$conditions <- "ExEvEs"

# Non-ExEvE days
non_exeves_prec_days <- exeves_all[is.na(event_80_95_id) & !is.na(prec),
                                    .N / (GRID_CELL_N * 0.5 * PERIOD_LENGTH),
                                    .(period, prec_day, month(date))]
setnames(non_exeves_prec_days, "V1", "value")
non_exeves_prec_days$conditions <- "Non-ExEvEs"

prec_days <- rbind(all_prec_days, exeves_prec_days, non_exeves_prec_days)
prec_days <- prec_days[, c(1, 3:2, 5:4)]

levels(prec_days$prec_day) <- c("Wet", "Dry")
levels(prec_days$period) <- c("Up to 2001", "After 2001")
setnames(prec_days, "prec_day", "Day class")

gg_wet_ratio <- ggplot(prec_days[conditions == 'ExEvEs']) +
  geom_col(aes(x = period, y = value, fill = `Day class`),
           position = position_fill(), width = 0.7) +
  facet_wrap(~ month) +
  xlab("Period") +
  ylab("Ratio") +
  scale_fill_manual(values = c('#a9cce0', '#fcc47c')) +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

#===============================================================================
# Combine & save
#===============================================================================
cat("Saving figure...\n")
ggarrange(gg_extreme_day, gg_mean_day, gg_wet_ratio,
          ncol = 1, labels = c("A", "B", "C"),
          legend = 'right', common.legend = TRUE)
ggsave(paste0(PATH_OUTPUT_FIGURES, "wet_days.png"), width = 9, height = 12)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "wet_days.png"), "\n")
rm(exeves, prec, exeves_prec, exeves_all); gc()
