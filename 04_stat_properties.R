# Statistical properties of extreme events – full comparison table
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/04_stat_properties.R

library(data.table)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
cat("Grid cells:", exeves[, uniqueN(grid_id)], "\n\n")

#===============================================================================
# OVERALL EVENT COMPARISON
#===============================================================================
event_comparison <- data.table(
  definition = factor(c("Q80/Q95", "Mean/Q95", "Mean/Q95*", "Q95", "Q95*", "Q80"),
                      levels = c("Q80/Q95", "Mean/Q95", "Mean/Q95*", "Q95", "Q95*", "Q80")))

event_comparison$frequency <- round(c(
  mean(exeves[!is.na(event_80_95_id), unique(event_80_95_id), .(grid_id)][, .N, grid_id]$N),
  mean(exeves[!is.na(event_id),       unique(event_id),       .(grid_id)][, .N, grid_id]$N),
  mean(exeves[!is.na(event_qr_id),    unique(event_qr_id),    .(grid_id)][, .N, grid_id]$N),
  mean(exeves[!is.na(extreme_id),     unique(extreme_id),     .(grid_id)][, .N, grid_id]$N),
  mean(exeves[!is.na(extreme_qr_id),  unique(extreme_qr_id),  .(grid_id)][, .N, grid_id]$N),
  mean(exeves[!is.na(event_80_id),    unique(event_80_id),    .(grid_id)][, .N, grid_id]$N)
) / PERIOD_LENGTH, 0)

event_comparison$duration_mean <- round(c(
  mean(exeves[!is.na(event_80_95_id), .N, .(grid_id, event_80_95_id)]$N),
  mean(exeves[!is.na(event_id),       .N, .(grid_id, event_id)]$N),
  mean(exeves[!is.na(event_qr_id),    .N, .(grid_id, event_qr_id)]$N),
  mean(exeves[!is.na(extreme_id),     .N, .(grid_id, extreme_id)]$N),
  mean(exeves[!is.na(extreme_qr_id),  .N, .(grid_id, extreme_qr_id)]$N),
  mean(exeves[!is.na(event_80_id),    .N, .(grid_id, event_80_id)]$N)
), 1)

event_comparison$duration_max <- c(
  max(exeves[!is.na(event_80_95_id), .N, .(grid_id, event_80_95_id)]$N),
  max(exeves[!is.na(event_id),       .N, .(grid_id, event_id)]$N),
  max(exeves[!is.na(event_qr_id),    .N, .(grid_id, event_qr_id)]$N),
  max(exeves[!is.na(extreme_id),     .N, .(grid_id, extreme_id)]$N),
  max(exeves[!is.na(extreme_qr_id),  .N, .(grid_id, extreme_qr_id)]$N),
  max(exeves[!is.na(event_80_id),    .N, .(grid_id, event_80_id)]$N))

event_comparison$extreme_n <- round(c(
  mean(exeves[!is.na(extreme_id),    .N, .(grid_id, event_80_95_id)]$N),
  mean(exeves[!is.na(extreme_id),    .N, .(grid_id, event_id)]$N),
  mean(exeves[!is.na(extreme_qr_id), .N, .(grid_id, event_qr_id)]$N),
  mean(exeves[!is.na(extreme_id),    .N, .(grid_id, extreme_id)]$N),
  mean(exeves[!is.na(extreme_qr_id), .N, .(grid_id, extreme_qr_id)]$N),
  mean(exeves[!is.na(extreme_id),    .N, .(grid_id, event_80_id)]$N)
), 1)

event_comparison$extreme_n_max <- round(c(
  max(exeves[!is.na(extreme_id),    .N, .(grid_id, event_80_95_id)]$N),
  max(exeves[!is.na(extreme_id),    .N, .(grid_id, event_id)]$N),
  max(exeves[!is.na(extreme_qr_id), .N, .(grid_id, event_qr_id)]$N),
  max(exeves[!is.na(extreme_id),    .N, .(grid_id, extreme_id)]$N),
  max(exeves[!is.na(extreme_qr_id), .N, .(grid_id, extreme_qr_id)]$N),
  max(exeves[!is.na(extreme_id),    .N, .(grid_id, event_80_id)]$N)
), 1)

event_comparison$intensity_mean <- round(c(
  mean(exeves[!is.na(event_80_95_id), mean(value), .(grid_id, event_80_95_id)]$V1),
  mean(exeves[!is.na(event_id),       mean(value), .(grid_id, event_id)]$V1),
  mean(exeves[!is.na(event_qr_id),    mean(value), .(grid_id, event_qr_id)]$V1),
  mean(exeves[!is.na(extreme_id),     mean(value), .(grid_id, extreme_id)]$V1),
  mean(exeves[!is.na(extreme_qr_id),  mean(value), .(grid_id, extreme_qr_id)]$V1),
  mean(exeves[!is.na(event_80_id),    mean(value), .(grid_id, event_80_id)]$V1)
), 1)

event_comparison$intensity_max <- round(c(
  max(exeves[!is.na(event_80_95_id), mean(value), .(grid_id, event_80_95_id)]$V1),
  max(exeves[!is.na(event_id),       mean(value), .(grid_id, event_id)]$V1),
  max(exeves[!is.na(event_qr_id),    mean(value), .(grid_id, event_qr_id)]$V1),
  max(exeves[!is.na(extreme_id),     mean(value), .(grid_id, extreme_id)]$V1),
  max(exeves[!is.na(extreme_qr_id),  mean(value), .(grid_id, extreme_qr_id)]$V1),
  max(exeves[!is.na(event_80_id),    mean(value), .(grid_id, event_80_id)]$V1)
), 1)

event_comparison$severity_mean <- round(c(
  mean(exeves[!is.na(event_80_95_id), sum(value), .(grid_id, event_80_95_id)]$V1),
  mean(exeves[!is.na(event_id),       sum(value), .(grid_id, event_id)]$V1),
  mean(exeves[!is.na(event_qr_id),    sum(value), .(grid_id, event_qr_id)]$V1),
  mean(exeves[!is.na(extreme_id),     sum(value), .(grid_id, extreme_id)]$V1),
  mean(exeves[!is.na(extreme_qr_id),  sum(value), .(grid_id, extreme_qr_id)]$V1),
  mean(exeves[!is.na(event_80_id),    sum(value), .(grid_id, event_80_id)]$V1)
), 1)

event_comparison$severity_max <- round(c(
  max(exeves[!is.na(event_80_95_id), sum(value), .(grid_id, event_80_95_id)]$V1),
  max(exeves[!is.na(event_id),       sum(value), .(grid_id, event_id)]$V1),
  max(exeves[!is.na(event_qr_id),    sum(value), .(grid_id, event_qr_id)]$V1),
  max(exeves[!is.na(extreme_id),     sum(value), .(grid_id, extreme_id)]$V1),
  max(exeves[!is.na(extreme_qr_id),  sum(value), .(grid_id, extreme_qr_id)]$V1),
  max(exeves[!is.na(event_80_id),    sum(value), .(grid_id, event_80_id)]$V1)
), 1)

write.csv(event_comparison, paste0(PATH_OUTPUT_TABLES, region, '_definition_properties.csv'),
          row.names = FALSE)
cat("Overall properties saved.\n")
print(event_comparison)

#===============================================================================
# PERIOD COMPARISON – helper
#===============================================================================
period_stats <- function(per) {
  dt <- data.table(
    definition = factor(c("Q80/Q95", "Mean/Q95", "Mean/Q95*", "Q95", "Q95*", "Q80"),
                        levels = c("Q80/Q95", "Mean/Q95", "Mean/Q95*", "Q95", "Q95*", "Q80")))
  dt$frequency <- round(c(
    mean(exeves[!is.na(event_80_95_id) & period == per, unique(event_80_95_id), .(grid_id)][, .N, grid_id]$N),
    mean(exeves[!is.na(event_id)       & period == per, unique(event_id),       .(grid_id)][, .N, grid_id]$N),
    mean(exeves[!is.na(event_qr_id)    & period == per, unique(event_qr_id),    .(grid_id)][, .N, grid_id]$N),
    mean(exeves[!is.na(extreme_id)     & period == per, unique(extreme_id),     .(grid_id)][, .N, grid_id]$N),
    mean(exeves[!is.na(extreme_qr_id)  & period == per, unique(extreme_qr_id),  .(grid_id)][, .N, grid_id]$N),
    mean(exeves[!is.na(event_80_id)    & period == per, unique(event_80_id),    .(grid_id)][, .N, grid_id]$N)
  ) / (0.5 * PERIOD_LENGTH), 0)

  dt$duration_mean <- round(c(
    mean(exeves[!is.na(event_80_95_id) & period == per, .N, .(grid_id, event_80_95_id)]$N),
    mean(exeves[!is.na(event_id)       & period == per, .N, .(grid_id, event_id)]$N),
    mean(exeves[!is.na(event_qr_id)    & period == per, .N, .(grid_id, event_qr_id)]$N),
    mean(exeves[!is.na(extreme_id)     & period == per, .N, .(grid_id, extreme_id)]$N),
    mean(exeves[!is.na(extreme_qr_id)  & period == per, .N, .(grid_id, extreme_qr_id)]$N),
    mean(exeves[!is.na(event_80_id)    & period == per, .N, .(grid_id, event_80_id)]$N)
  ), 1)

  dt$duration_max <- c(
    max(exeves[!is.na(event_80_95_id) & period == per, .N, .(grid_id, event_80_95_id)]$N),
    max(exeves[!is.na(event_id)       & period == per, .N, .(grid_id, event_id)]$N),
    max(exeves[!is.na(event_qr_id)    & period == per, .N, .(grid_id, event_qr_id)]$N),
    max(exeves[!is.na(extreme_id)     & period == per, .N, .(grid_id, extreme_id)]$N),
    max(exeves[!is.na(extreme_qr_id)  & period == per, .N, .(grid_id, extreme_qr_id)]$N),
    max(exeves[!is.na(event_80_id)    & period == per, .N, .(grid_id, event_80_id)]$N))

  dt$intensity_mean <- round(c(
    mean(exeves[!is.na(event_80_95_id) & period == per, mean(value), .(grid_id, event_80_95_id)]$V1),
    mean(exeves[!is.na(event_id)       & period == per, mean(value), .(grid_id, event_id)]$V1),
    mean(exeves[!is.na(event_qr_id)    & period == per, mean(value), .(grid_id, event_qr_id)]$V1),
    mean(exeves[!is.na(extreme_id)     & period == per, mean(value), .(grid_id, extreme_id)]$V1),
    mean(exeves[!is.na(extreme_qr_id)  & period == per, mean(value), .(grid_id, extreme_qr_id)]$V1),
    mean(exeves[!is.na(event_80_id)    & period == per, mean(value), .(grid_id, event_80_id)]$V1)
  ), 1)

  dt$severity_mean <- round(c(
    mean(exeves[!is.na(event_80_95_id) & period == per, sum(value), .(grid_id, event_80_95_id)]$V1),
    mean(exeves[!is.na(event_id)       & period == per, sum(value), .(grid_id, event_id)]$V1),
    mean(exeves[!is.na(event_qr_id)    & period == per, sum(value), .(grid_id, event_qr_id)]$V1),
    mean(exeves[!is.na(extreme_id)     & period == per, sum(value), .(grid_id, extreme_id)]$V1),
    mean(exeves[!is.na(extreme_qr_id)  & period == per, sum(value), .(grid_id, extreme_qr_id)]$V1),
    mean(exeves[!is.na(event_80_id)    & period == per, sum(value), .(grid_id, event_80_id)]$V1)
  ), 1)

  return(dt)
}

event_comparison_1 <- period_stats('up_to_2001')
event_comparison_2 <- period_stats('after_2001')

# Numeric columns only for ratio
num_cols <- names(event_comparison_1)[sapply(event_comparison_1, is.numeric)]
event_comparison_change <- copy(event_comparison_2)
event_comparison_change[, (num_cols) := lapply(.SD, function(x) round(x / event_comparison_1[[cur_column()]], 2)),
                        .SDcols = num_cols]

write.csv(event_comparison_change, paste0(PATH_OUTPUT_TABLES, region, '_definition_change.csv'),
          row.names = FALSE)

cat("\nPeriod change ratios saved.\n")
print(event_comparison_change)

rm(exeves); gc()
cat("\nStatistical properties analysis complete!\n")
