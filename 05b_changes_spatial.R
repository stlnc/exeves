# Spatial change analysis: severity, intensity, frequency, duration between periods
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/05b_changes_spatial.R
# Note: radiation data not available; using evaporation + precipitation only

library(data.table)
library(ggplot2)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

cat("Loading data...\n")
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))
exeves    <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec      <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))

# Merge precipitation – keyed join (no full copy)
setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.value, on = .(grid_id, date)]

# Add spatial coordinates via keyed update
setkey(evap_grid, grid_id)
exeves[evap_grid, `:=`(lon = i.lon, lat = i.lat), on = "grid_id"]
rm(evap_grid, prec); gc()

#===============================================================================
# SEVERITY: All values
#===============================================================================
cat("Computing severity changes...\n")

evap_severity_period <- exeves[, .(value = sum(value)), .(lon, lat, period)]
evap_severity_period[, diff_value := diff(value), by = .(lon, lat)]
evap_severity_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
evap_severity_period$variable <- "E: All days"

prec_severity_period <- exeves[, .(value = sum(prec, na.rm = TRUE)), .(lon, lat, period)]
prec_severity_period[, diff_value := diff(value), by = .(lon, lat)]
prec_severity_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
prec_severity_period$variable <- "P: All days"

#===============================================================================
# SEVERITY: ExEvEs
#===============================================================================
event_evap_severity_period <- exeves[!is.na(event_id), .(value = sum(value)), .(lon, lat, period)]
event_evap_severity_period <- event_evap_severity_period[event_evap_severity_period[, .N, by = .(lon, lat)][N == 2], on = .(lon, lat)]
event_evap_severity_period[, diff_value := diff(value), by = .(lon, lat)]
event_evap_severity_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
event_evap_severity_period$variable <- "Evaporation (ExEvEs)"

event_prec_severity_period <- exeves[!is.na(event_id), .(value = sum(prec, na.rm = TRUE)), .(lon, lat, period)]
event_prec_severity_period <- event_prec_severity_period[event_prec_severity_period[, .N, by = .(lon, lat)][N == 2], on = .(lon, lat)]
event_prec_severity_period[, diff_value := diff(value), by = .(lon, lat)]
event_prec_severity_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
event_prec_severity_period$variable <- "Precipitation (ExEvEs)"

#===============================================================================
# INTENSITY
#===============================================================================
cat("Computing intensity changes...\n")
event_evap_intensity_period <- exeves[!is.na(event_id), .(value = round(mean(value), 2)), .(lon, lat, period)]
event_evap_intensity_period <- event_evap_intensity_period[event_evap_intensity_period[, .N, by = .(lon, lat)][N == 2], on = .(lon, lat)]
event_evap_intensity_period[, diff_value := diff(value), by = .(lon, lat)]
event_evap_intensity_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
event_evap_intensity_period$variable <- "Intensity (E)"

event_prec_intensity_period <- exeves[!is.na(event_id), .(value = round(mean(prec, na.rm = TRUE), 2)), .(lon, lat, period)]
event_prec_intensity_period <- event_prec_intensity_period[event_prec_intensity_period[, .N, by = .(lon, lat)][N == 2], on = .(lon, lat)]
event_prec_intensity_period[, diff_value := diff(value), by = .(lon, lat)]
event_prec_intensity_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
event_prec_intensity_period$variable <- "Intensity (P)"

#===============================================================================
# FREQUENCY
#===============================================================================
cat("Computing frequency changes...\n")
event_frequency_period <- exeves[!is.na(event_id), .(value = .N), .(lon, lat, period)]
event_frequency_period <- event_frequency_period[event_frequency_period[, .N, by = .(lon, lat)][N == 2], on = .(lon, lat)]
event_frequency_period[, diff_value := diff(value), by = .(lon, lat)]
event_frequency_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
event_frequency_period$variable <- factor("ExEvEs Frequency")

#===============================================================================
# DURATION
#===============================================================================
cat("Computing duration changes...\n")
event_duration_period <- exeves[!is.na(event_id), .(value = .N), .(event_id, lon, lat, period)]
event_duration_period <- event_duration_period[, .(value = mean(value)), .(lon, lat, period)]
event_duration_period <- event_duration_period[event_duration_period[, .N, by = .(lon, lat)][N == 2], on = .(lon, lat)]
event_duration_period[, diff_value := diff(value), by = .(lon, lat)]
event_duration_period[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
event_duration_period$variable <- factor("ExEvEs Duration")

#===============================================================================
# COMBINE AND SAVE
#===============================================================================
exeve_properties_change <- rbind(
  event_evap_severity_period, event_prec_severity_period,
  event_evap_intensity_period, event_prec_intensity_period,
  event_frequency_period, event_duration_period)
saveRDS(exeve_properties_change, file = paste0(PATH_OUTPUT, 'spatial_changes.rds'))

cat("Spatial changes saved.\n")

# Quick validation plot
ggplot(event_evap_severity_period[period == "up_to_2001"]) +
  geom_tile(aes(lon, lat, fill = ratio)) +
  scale_fill_gradient2(low = "navyblue", mid = "grey90", high = "darkred", midpoint = 1) +
  theme_minimal()
ggsave(paste0(PATH_OUTPUT_FIGURES, "spatial_changes_validation.png"), width = 10, height = 8)

rm(exeves); gc()
cat("Spatial change analysis complete.\n")
