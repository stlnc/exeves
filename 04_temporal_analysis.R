# Temporal trends and seasonal patterns
# Analyzes how extreme events change over time

library(data.table)
library(ggplot2)
library(lubridate)

# Memory optimization settings
setDTthreads(1)  # Reduce parallel threads to save memory

# Load paths and data
load("paths.Rdata")
region <- 'europe_med'

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))

# Annual time series
cat("Calculating annual statistics...\n")
annual_stats <- exeves[!is.na(event_80_95_id), .(
  n_events = uniqueN(event_80_95_id),
  total_days = as.numeric(.N),
  mean_duration = mean(as.numeric(.N)),
  total_evap = sum(value)
), .(year = year(date), grid_id)][, .(
  events_per_cell = mean(n_events),
  total_event_days = sum(total_days),
  mean_duration = mean(mean_duration),
  total_evap = sum(total_evap)
), year]
gc()

# Monthly climatology
monthly_stats <- exeves[!is.na(event_80_95_id), .(
  n_events = uniqueN(event_80_95_id)
), .(month = month(date, label = TRUE), grid_id)][, .(
  events_per_cell = mean(n_events)
), month]
gc()

# Seasonal statistics
seasonal_stats <- exeves[!is.na(event_80_95_id), .(
  n_events = uniqueN(event_80_95_id)
), .(season, grid_id, period)][, .(
  events_per_cell = mean(n_events),
  total_events = sum(n_events)
), .(season, period)]
gc()

# Create plots
cat("Creating plots...\n")

# Plot 1: Annual time series
p1 <- ggplot(annual_stats, aes(x = year, y = events_per_cell)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = 2001.5, linetype = "dotted", color = "gray40") +
  labs(title = "Annual Extreme Event Frequency",
       subtitle = "Mean events per grid cell (Q80/Q95 definition)",
       x = "Year", y = "Events per cell") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "timeseries_annual.png"), 
       p1, width = 10, height = 6, dpi = 300)
rm(p1); gc()

# Plot 2: Monthly climatology
p2 <- ggplot(monthly_stats, aes(x = month, y = events_per_cell)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(title = "Monthly Distribution of Extreme Events",
       subtitle = "Mean events per grid cell by month (1981-2022)",
       x = "Month", y = "Events per cell") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "climatology_monthly.png"), 
       p2, width = 10, height = 6, dpi = 300)
rm(p2); gc()

# Plot 3: Seasonal comparison between periods
seasonal_comparison <- dcast(seasonal_stats, season ~ period, value.var = "events_per_cell")
seasonal_comparison_long <- melt(seasonal_comparison, id.vars = "season", 
                                 variable.name = "period", value.name = "events")

p3 <- ggplot(seasonal_comparison_long, aes(x = season, y = events, fill = period)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("up_to_2001" = "steelblue", "after_2001" = "darkorange2"),
                    labels = c("1981-2001", "2002-2022")) +
  labs(title = "Seasonal Event Frequency by Period",
       subtitle = "Mean events per grid cell",
       x = "Season", y = "Events per cell", fill = "Period") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "seasonal_comparison.png"), 
       p3, width = 10, height = 6, dpi = 300)
rm(p3, seasonal_comparison, seasonal_comparison_long); gc()

# Plot 4: Event duration over time
duration_annual <- exeves[!is.na(event_80_95_id), .N, 
                         .(event_80_95_id, grid_id, year = year(date))][, .(
  mean_duration = mean(N),
  median_duration = median(N),
  q25 = quantile(N, 0.25),
  q75 = quantile(N, 0.75)
), year]

p4 <- ggplot(duration_annual, aes(x = year, y = mean_duration)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = 2001.5, linetype = "dotted", color = "gray40") +
  labs(title = "Event Duration Trends",
       subtitle = "Mean duration (solid line) with IQR (shaded area)",
       x = "Year", y = "Duration (days)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(PATH_OUTPUT_FIGURES, "timeseries_duration.png"), 
       p4, width = 10, height = 6, dpi = 300)
rm(p4, duration_annual); gc()

# Statistical tests
cat("\n=== Temporal Trend Analysis ===\n")

# Linear trend in annual frequency
trend_model <- lm(events_per_cell ~ year, data = annual_stats)
cat("\nAnnual frequency trend:\n")
cat("  Slope:", round(coef(trend_model)[2], 4), "events/cell/year\n")
cat("  p-value:", format.pval(summary(trend_model)$coefficients[2,4], digits = 3), "\n")
cat("  R-squared:", round(summary(trend_model)$r.squared, 3), "\n")

# Period comparison
period_comparison <- exeves[!is.na(event_80_95_id), uniqueN(event_80_95_id), 
                            .(grid_id, period)][, .(mean_events = mean(V1)), period]
cat("\nMean events per cell by period:\n")
print(period_comparison)

period_increase <- period_comparison$mean_events[2] / period_comparison$mean_events[1]
cat("Period increase factor:", round(period_increase, 3), "\n")

# Seasonal summary
cat("\n=== Seasonal Analysis ===\n")
print(seasonal_stats)

# Save summary table
summary_table <- data.table(
  Metric = c("Total years", "Mean events/cell/year", "Trend (events/cell/year)", 
             "Period 1 mean", "Period 2 mean", "Change factor"),
  Value = c(
    annual_stats[, .N],
    round(mean(annual_stats$events_per_cell), 2),
    round(coef(trend_model)[2], 4),
    round(period_comparison$mean_events[1], 2),
    round(period_comparison$mean_events[2], 2),
    round(period_increase, 3)
  )
)

write.csv(summary_table, paste0(PATH_OUTPUT_TABLES, "temporal_summary.csv"), row.names = FALSE)

cat("\nPlots saved to:", PATH_OUTPUT_FIGURES, "\n")

rm(exeves, annual_stats, monthly_stats, seasonal_stats, trend_model, period_comparison, summary_table)
gc()

cat("Temporal analysis complete!\n")
