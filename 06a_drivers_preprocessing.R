# Drivers preprocessing: merge ExEvE data with precipitation
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/06a_drivers_preprocessing.R
# Note: radiation/temperature/sensible heat data not available for this region;
#       precipitation is used as the available co-variable (driver)

library(data.table)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec   <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))

# Merge precipitation – keyed update-join (avoids a full copy)
setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.value, on = .(grid_id, date)]
exeves_drivers <- exeves
rm(exeves, prec); gc()

# Classify conditions
exeves_drivers[, conditions := ordered('ExEvE')]
exeves_drivers[is.na(event_80_95_id), conditions := ordered('non-ExEvE')]

# Event day and duration (within each event)
exeves_drivers[conditions == 'ExEvE',
               event_day := seq_len(.N), by = .(event_80_95_id, grid_id)]
exeves_drivers[conditions == 'ExEvE',
               event_duration := .N, by = .(event_80_95_id, grid_id)]

# Save full version
saveRDS(exeves_drivers, paste0(PATH_OUTPUT_DATA, region, '_exeves_drivers_all.rds'))

# Save trimmed version
exeves_drivers <- exeves_drivers[, .(grid_id, date, conditions, event_day, event_duration,
                                      evap = value, std_value, prec)]
saveRDS(exeves_drivers, paste0(PATH_OUTPUT_DATA, region, '_exeves_drivers.rds'))

cat("Drivers preprocessing complete.\n")
cat("  Saved:", paste0(region, '_exeves_drivers.rds'), "\n")
rm(exeves_drivers); gc()
