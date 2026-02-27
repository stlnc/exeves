# Spatial visualization of extreme events
# Creates maps showing event frequency, duration, and severity

library(data.table)
library(ggplot2)
library(sf)

# Load paths and data
load("paths.Rdata")
region <- 'europe_med'

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

# Calculate spatial statistics for Q80/Q95 definition
cat("Calculating spatial statistics...\n")
spatial_stats <- exeves[!is.na(event_80_95_id), .(
  n_events = uniqueN(event_80_95_id),
  mean_duration = mean(.N),
  total_evap = sum(value),
  mean_intensity = mean(value)
), .(grid_id, event_80_95_id)][, .(
  event_frequency = .N,
  mean_duration = mean(mean_duration),
  mean_severity = mean(total_evap),
  mean_intensity = mean(mean_intensity)
), grid_id]

# Merge with coordinates
spatial_data <- merge(evap_grid, spatial_stats, by = "grid_id", all.x = TRUE)
spatial_data[is.na(event_frequency), event_frequency := 0]

# Period-specific statistics
spatial_period1 <- exeves[!is.na(event_80_95_id) & period == "up_to_2001", 
                          uniqueN(event_80_95_id), .(grid_id)]
setnames(spatial_period1, "V1", "events_period1")

spatial_period2 <- exeves[!is.na(event_80_95_id) & period == "after_2001", 
                          uniqueN(event_80_95_id), .(grid_id)]
setnames(spatial_period2, "V1", "events_period2")

spatial_change <- merge(spatial_period1, spatial_period2, by = "grid_id", all = TRUE)
spatial_change[is.na(events_period1), events_period1 := 0]
spatial_change[is.na(events_period2), events_period2 := 0]
spatial_change[, change := events_period2 - events_period1]
spatial_change[, ratio := ifelse(events_period1 > 0, events_period2 / events_period1, NA)]

spatial_change <- merge(evap_grid, spatial_change, by = "grid_id", all.x = TRUE)

cat("Creating plots...\n")

# Plot 1: Event frequency
p1 <- ggplot(spatial_data, aes(x = lon, y = lat, fill = event_frequency)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Events") +
  coord_quickmap() +
  labs(title = "Extreme Event Frequency (1981-2022)",
       subtitle = "Q80/Q95 definition",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_line(color = "gray90")
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "map_event_frequency.png"), 
       p1, width = 10, height = 8, dpi = 300)

# Plot 2: Mean event duration
p2 <- ggplot(spatial_data[event_frequency > 0], aes(x = lon, y = lat, fill = mean_duration)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "Days") +
  coord_quickmap() +
  labs(title = "Mean Event Duration",
       subtitle = "Average days per extreme event",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_line(color = "gray90")
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "map_event_duration.png"), 
       p2, width = 10, height = 8, dpi = 300)

# Plot 3: Change between periods
p3 <- ggplot(spatial_change, aes(x = lon, y = lat, fill = change)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Change"
  ) +
  coord_quickmap() +
  labs(title = "Change in Event Frequency",
       subtitle = "Difference: 2002-2022 minus 1981-2001",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_line(color = "gray90")
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "map_event_change.png"), 
       p3, width = 10, height = 8, dpi = 300)

# Plot 4: Ratio of change
p4 <- ggplot(spatial_change[!is.na(ratio) & events_period1 >= 5], 
             aes(x = lon, y = lat, fill = ratio)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "dodgerblue", mid = "gray90", high = "darkorange3", midpoint = 1,
    name = "Ratio",
    limits = c(0.5, 1.5)
  ) +
  coord_quickmap() +
  labs(title = "Ratio of Event Frequency Change",
       subtitle = "2002-2022 / 1981-2001 (cells with ≥5 events in first period)",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_line(color = "gray90")
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "map_event_ratio.png"), 
       p4, width = 10, height = 8, dpi = 300)

# Summary statistics
cat("\n=== Spatial Summary ===\n")
cat("Total grid cells:", spatial_data[, .N], "\n")
cat("Cells with events:", spatial_data[event_frequency > 0, .N], "\n")
cat("Mean events per cell:", round(mean(spatial_data$event_frequency, na.rm = TRUE), 1), "\n")
cat("\nFrequency distribution:\n")
print(summary(spatial_data$event_frequency))

cat("\n=== Change Summary ===\n")
cat("Cells with increase:", spatial_change[change > 0, .N], "\n")
cat("Cells with decrease:", spatial_change[change < 0, .N], "\n")
cat("Cells with no change:", spatial_change[change == 0, .N], "\n")
cat("\nMean change:", round(mean(spatial_change$change, na.rm = TRUE), 2), "events\n")

cat("\nMaps saved to:", PATH_OUTPUT_FIGURES, "\n")
cat("Visualization complete!\n")
