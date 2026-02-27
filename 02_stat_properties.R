# Statistical properties of extreme events
# Compares different ExEvE definitions

library(data.table)

# Load paths and data
load("paths.Rdata")
region <- 'europe_med'

cat("Loading exeves data for", region, "...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))

cat("Grid cells:", exeves[, uniqueN(grid_id)], "\n")
cat("Date range:", as.character(min(exeves$date)), "to", as.character(max(exeves$date)), "\n\n")

# Compare event definitions
event_comparison <- data.table(
  definition = factor(c("Q80/Q95", "Mean/Q95", "Mean/Q95*", "Q80"), 
                     levels = c("Q80/Q95", "Mean/Q95", "Mean/Q95*", "Q80"))
)

# Event frequency (mean events per grid cell)
event_comparison$frequency <- round(c(
  mean(exeves[!is.na(event_80_95_id), uniqueN(event_80_95_id), grid_id]$V1),
  mean(exeves[!is.na(event_id), uniqueN(event_id), grid_id]$V1),
  mean(exeves[!is.na(event_qr_id), uniqueN(event_qr_id), grid_id]$V1),
  mean(exeves[!is.na(event_80_id), uniqueN(event_80_id), grid_id]$V1)
), 1)

# Event duration (mean days per event)
event_comparison$duration <- round(c(
  mean(exeves[!is.na(event_80_95_id), .N, .(event_80_95_id, grid_id)]$N),
  mean(exeves[!is.na(event_id), .N, .(event_id, grid_id)]$N),
  mean(exeves[!is.na(event_qr_id), .N, .(event_qr_id, grid_id)]$N),
  mean(exeves[!is.na(event_80_id), .N, .(event_80_id, grid_id)]$N)
), 1)

# Maximum event duration
event_comparison$max_duration <- c(
  max(exeves[!is.na(event_80_95_id), .N, .(event_80_95_id, grid_id)]$N),
  max(exeves[!is.na(event_id), .N, .(event_id, grid_id)]$N),
  max(exeves[!is.na(event_qr_id), .N, .(event_qr_id, grid_id)]$N),
  max(exeves[!is.na(event_80_id), .N, .(event_80_id, grid_id)]$N)
)

# Event severity (mean evaporation sum per event, mm)
event_comparison$severity <- round(c(
  mean(exeves[!is.na(event_80_95_id), sum(value), .(event_80_95_id, grid_id)]$V1),
  mean(exeves[!is.na(event_id), sum(value), .(event_id, grid_id)]$V1),
  mean(exeves[!is.na(event_qr_id), sum(value), .(event_qr_id, grid_id)]$V1),
  mean(exeves[!is.na(event_80_id), sum(value), .(event_80_id, grid_id)]$V1)
), 1)

# Event intensity (mean evaporation rate during events, mm/day)
event_comparison$intensity <- round(c(
  mean(exeves[!is.na(event_80_95_id), mean(value), .(event_80_95_id, grid_id)]$V1),
  mean(exeves[!is.na(event_id), mean(value), .(event_id, grid_id)]$V1),
  mean(exeves[!is.na(event_qr_id), mean(value), .(event_qr_id, grid_id)]$V1),
  mean(exeves[!is.na(event_80_id), mean(value), .(event_80_id, grid_id)]$V1)
), 2)

print(event_comparison)

# Save results
write.csv(event_comparison, paste0(PATH_OUTPUT_TABLES, region, '_event_properties.csv'), row.names = FALSE)
cat("\nResults saved to:", paste0(PATH_OUTPUT_TABLES, region, '_event_properties.csv\n'))

# Period comparison (1981-2001 vs 2002-2022)
cat("\n=== Temporal Changes ===\n")
cat("\nMean events per grid cell by period:\n")
period_freq <- exeves[!is.na(event_80_95_id), uniqueN(event_80_95_id), .(grid_id, period)][
  , .(mean_events = mean(V1)), period]
print(period_freq)

cat("\nMean event duration by period (days):\n")
period_dur <- exeves[!is.na(event_80_95_id), .N, .(event_80_95_id, grid_id, period)][
  , .(mean_duration = mean(N)), period]
print(period_dur)

# Seasonal distribution
cat("\n=== Seasonal Distribution ===\n")
cat("\nEvents by season (Q80/Q95 definition):\n")
seasonal_events <- exeves[!is.na(event_80_95_id), uniqueN(event_80_95_id), .(grid_id, season)][
  , .(n_events = sum(V1)), season]
seasonal_events[, percentage := round(100 * n_events / sum(n_events), 1)]
print(seasonal_events)

# Spatial summary
cat("\n=== Spatial Summary ===\n")
spatial_summary <- exeves[!is.na(event_80_95_id), .(
  n_events = uniqueN(event_80_95_id),
  mean_duration = mean(.N), 
  total_severity = sum(value)
), grid_id]

cat("Grid cells with events:", spatial_summary[, .N], "\n")
cat("Mean events per cell:", round(mean(spatial_summary$n_events), 1), "\n")
cat("Range of events per cell:", min(spatial_summary$n_events), "to", max(spatial_summary$n_events), "\n")

# Example events (longest events)
cat("\n=== Example Events (longest 5 events) ===\n")
longest_events <- exeves[!is.na(event_80_95_id), .(
  duration = .N,
  severity = sum(value),
  intensity = mean(value),
  start_date = min(date),
  end_date = max(date)
), .(event_80_95_id, grid_id)][order(-duration)][1:5]

print(longest_events)

cat("\nStatistical analysis complete!\n")
