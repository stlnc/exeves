# Preprocessing evaporation data for Europe/Mediterranean
# Identifies extreme evaporation events (ExEvEs)

# Load required packages
library(data.table)
library(ncdf4)
library(lubridate)

# Memory optimization settings
setDTthreads(1)  # Reduce parallel threads to save memory

# Load paths
load("paths.Rdata")

# Configuration
region <- 'europe_med'

# Time periods
START_PERIOD_1 <- as.Date("1981-1-1") 
END_PERIOD_1 <- as.Date("2001-12-31")
END_PERIOD_2 <- as.Date("2022-12-31")
PERIOD_LENGTH <- round(as.numeric((END_PERIOD_2 - START_PERIOD_1) / 365.25), 0)

# Thresholds
EXTREMES_THRES <- 0.95  # 95th percentile for extremes
LOW_THRES <- 0.80       # 80th percentile for event start
SUB_PERIOD_YEARS <- 0.5 * PERIOD_LENGTH
DAYS_IN_YEAR <- 365

cat("Processing", region, "data...\n")
cat("Period:", as.character(START_PERIOD_1), "to", as.character(END_PERIOD_2), "\n")

# Read NetCDF file
nc_file <- paste0("gleam_e_mm_", region, ".nc")
cat("Reading:", nc_file, "\n")

nc <- nc_open(nc_file)
evap_array <- ncvar_get(nc, "E")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")

# Get time attributes to convert to dates
time_units <- ncatt_get(nc, "time", "units")$value
# Typical format: "days since 1980-01-01" or "hours since..."
time_start <- as.Date(sub(".*since ", "", time_units))

if (grepl("hours", time_units)) {
  dates <- time_start + time / 24
} else {
  dates <- time_start + time
}

nc_close(nc)
rm(nc)  # Remove nc object
gc()  # Force garbage collection

cat("Dimensions: lon =", length(lon), ", lat =", length(lat), ", time =", length(dates), "\n")

# Convert to data.table
cat("Converting to data.table format...\n")
evap_list <- list()

for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    evap_ts <- evap_array[i, j, ]
    if (!all(is.na(evap_ts))) {  # Skip all-NA pixels
      evap_list[[length(evap_list) + 1]] <- data.table(
        lon = lon[i],
        lat = lat[j],
        date = dates,
        value = evap_ts
      )
    }
  }
}

rm(evap_array)  # Remove large array after loop
gc()  # Force garbage collection

evap <- rbindlist(evap_list)
rm(evap_list)  # Remove intermediate list
gc()  # Force garbage collection
evap <- evap[!is.na(value)]  # Remove NA values
evap <- evap[order(lon, lat, date)]  # Order data

cat("Total observations:", nrow(evap), "\n")

# Create grid IDs
evap[, grid_id := .GRP, by = list(lat, lon)]
evap_grid <- unique(evap[, .(lon, lat, grid_id)])
cat("Number of grid cells:", evap_grid[, .N], "\n")

# Save grid
saveRDS(evap_grid, paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

# Keep only grid_id, date, value
evap <- evap[, .(grid_id, date, value)]

# Save raw evap data
cat("Saving evaporation grid data...\n")
saveRDS(evap, paste0(PATH_OUTPUT_DATA, region, '_evap_grid.rds'))
gc()  # Force garbage collection

## Create Pentads (5-day periods)
cat("\nCreating pentad statistics...\n")
pentads <- evap[, .(grid_id, date, value)]
pentads[, pentad := ceiling((yday(date) - leap_year(year(date)) * (yday(date) > 59)) / 5)]
pentads[, std_value := (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE), by = .(pentad, grid_id)]
pentads[, pentad_std_q95 := quantile(std_value, EXTREMES_THRES, na.rm = TRUE), by = grid_id]
pentads[, pentad_std_q80 := quantile(std_value, LOW_THRES, na.rm = TRUE), by = grid_id]
gc()  # Force garbage collection after large computation

# Additional quantile definitions for alternative thresholds
cat("Computing additional quantile thresholds...\n")
pentads[, pentad_median := quantile(std_value, 0.50, na.rm = TRUE), by = grid_id]
pentads[, pentad_std_q95_alt := quantile(std_value, EXTREMES_THRES, na.rm = TRUE), by = grid_id]
gc()  # Force garbage collection

## Identify Extreme Events
cat("\nIdentifying extreme events (memory-optimized with chunking)...\n")

# Process in chunks if dataset is very large (>50M rows)
n_obs <- nrow(evap)
CHUNK_SIZE_ROWS <- 20000000  # Process 20M rows at a time
process_in_chunks <- n_obs > 50000000

if (process_in_chunks) {
  cat(sprintf("Large dataset detected (%d rows). Processing in chunks to reduce memory...\n", n_obs))
  n_grids <- uniqueN(evap$grid_id)
  grids_per_chunk <- max(100, ceiling(n_grids * CHUNK_SIZE_ROWS / n_obs))
  cat(sprintf("Processing %d grids in chunks of ~%d grids each\n", n_grids, grids_per_chunk))
}

### Main definition: Mean/Q95
exeves <- merge(evap, pentads[, .(grid_id, date, std_value, pentad_std_q95, pentad_std_q80, pentad_median, pentad_std_q95_alt)], 
                all.x = TRUE, by = c("grid_id", "date"))
gc()  # Force garbage collection after large merge

cat("Computing event identifiers for Q80/Q95 definition...\n")
exeves[, evap_event := FALSE]
exeves[, value_above_low_thres := FALSE]
exeves[, extreme := FALSE]

# Event = period above low threshold (80th percentile) containing at least one extreme (95th percentile)
exeves[std_value > pentad_std_q80, value_above_low_thres := TRUE]
exeves[std_value > pentad_std_q95, extreme := TRUE]
setkeyv(exeves, c("grid_id", "date"))  # Set key for efficient operations
exeves[, above_low_thres_id := rleid(value_above_low_thres), by = grid_id]
exeves[, extreme_id := rleid(extreme), by = grid_id]

# Mark entire above-threshold period as event if it contains an extreme - memory-efficient version
exeves[extreme == TRUE, evap_event := TRUE, .(grid_id, above_low_thres_id)] 
periods_with_extreme <- unique(exeves[extreme == TRUE, .(grid_id, above_low_thres_id)])
setkey(exeves, grid_id, above_low_thres_id)
setkey(periods_with_extreme, grid_id, above_low_thres_id)
exeves[periods_with_extreme, evap_event := TRUE]
rm(periods_with_extreme)
exeves[, event_id := rleid(evap_event), .(grid_id)]
exeves[evap_event != TRUE, event_id := NA]
exeves[extreme != TRUE, extreme_id := NA]
gc()

# Add time periods
exeves[, period := ordered('up_to_2001')]
exeves[date > END_PERIOD_1, period := ordered('after_2001')]

# Add seasons
exeves[month(date) < 4, season := ordered("JFM")]
exeves[month(date) >= 4 & month(date) < 7, season := ordered("AMJ")]
exeves[month(date) >= 7 & month(date) < 10, season := ordered("JAS")]
exeves[month(date) >= 10, season := ordered("OND")]
gc()  # Force garbage collection

cat("Computing alternative definitions in-place...\n")

### Alternative definition: Median/Q95 - computed in-place
exeves[, evap_event_qr := FALSE]
exeves[, value_above_median := FALSE]
exeves[, extreme_qr := FALSE]
exeves[std_value > pentad_median, value_above_median := TRUE]
exeves[std_value > pentad_std_q95_alt, extreme_qr := TRUE]
exeves[, above_median_id := rleid(value_above_median), by = grid_id]

periods_with_extreme_qr <- unique(exeves[extreme_qr == TRUE, .(grid_id, above_median_id)])
setkey(exeves, grid_id, above_median_id)
setkey(periods_with_extreme_qr, grid_id, above_median_id)
exeves[periods_with_extreme_qr, evap_event_qr := TRUE]
rm(periods_with_extreme_qr)
exeves[, event_qr_id := rleid(evap_event_qr), .(grid_id)]
exeves[evap_event_qr != TRUE, event_qr_id := NA]
exeves[, extreme_qr_id := rleid(extreme_qr), .(grid_id)]
exeves[extreme_qr != TRUE, extreme_qr_id := NA]

# Clean up temporary columns
exeves[, c("evap_event_qr", "value_above_median", "extreme_qr", "above_median_id") := NULL]
gc()

### Q80/Q95 definition - computed in-place
exeves[, evap_event_80_95 := FALSE]
exeves[, value_above_q80 := FALSE]
exeves[, extreme_q95 := FALSE]
exeves[std_value > pentad_std_q80, value_above_q80 := TRUE]
exeves[std_value > pentad_std_q95, extreme_q95 := TRUE]
exeves[, above_q80_id := rleid(value_above_q80), by = grid_id]

periods_with_extreme_95 <- unique(exeves[extreme_q95 == TRUE, .(grid_id, above_q80_id)])
setkey(exeves, grid_id, above_q80_id)
setkey(periods_with_extreme_95, grid_id, above_q80_id)
exeves[periods_with_extreme_95, evap_event_80_95 := TRUE]
rm(periods_with_extreme_95)
exeves[, event_80_95_id := rleid(evap_event_80_95), .(grid_id)]
exeves[evap_event_80_95 != TRUE, event_80_95_id := NA]
exeves[, c("evap_event_80_95", "value_above_q80", "extreme_q95", "above_q80_id") := NULL]
gc()

### Simple Q80 definition - computed in-place
exeves[, evap_event_80 := FALSE]
exeves[std_value > pentad_std_q80, evap_event_80 := TRUE]
exeves[, event_80_id := rleid(evap_event_80), .(grid_id)]
exeves[evap_event_80 != TRUE, event_80_id := NA]

# Clean up all temporary columns
exeves[, c("evap_event", "value_above_low_thres", "extreme", "above_low_thres_id", "evap_event_80",
           "pentad_std_q95", "pentad_std_q80", "pentad_median", "pentad_std_q95_alt") := NULL]
gc()

# Keep only necessary columns and merge with coordinates
exeves <- exeves[, .(grid_id, date, season, period, value, std_value, 
                     event_id, extreme_id, event_qr_id, extreme_qr_id, 
                     event_80_95_id, event_80_id)]
                     
exeves <- merge(evap_grid, exeves, by = 'grid_id')
gc()

# Save results
cat("\nSaving results...\n")
saveRDS(pentads, paste0(PATH_OUTPUT_DATA, 'pentads_std_', region, '.rds'))
saveRDS(exeves, paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
gc()  # Force garbage collection after saving

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Total grid cells:", exeves[, uniqueN(grid_id)], "\n")
cat("Date range:", min(exeves$date), "to", max(exeves$date), "\n")
cat("\nEvent counts per grid cell (mean):\n")
cat("  Q80/Q95:", mean(exeves[!is.na(event_80_95_id), uniqueN(event_80_95_id), grid_id]$V1), "\n")
cat("  Mean/Q95:", mean(exeves[!is.na(event_id), uniqueN(event_id), grid_id]$V1), "\n")
cat("  Mean/Q95* (QR):", mean(exeves[!is.na(event_qr_id), uniqueN(event_qr_id), grid_id]$V1), "\n")
cat("  Q80 only:", mean(exeves[!is.na(event_80_id), uniqueN(event_80_id), grid_id]$V1), "\n")

cat("\nPreprocessing complete!\n")
rm(evap, pentads, exeves); gc()
