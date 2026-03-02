# Preprocessing: identify extreme evaporation events (ExEvEs) using multiple definitions
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/02_preprocessing.R

library(data.table)
library(lubridate)
library(quantreg)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

cat("Processing region:", region, "\n")

#===============================================================================
# 1. LOAD EVAPORATION DATA
#===============================================================================
evap <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_evap_grid.rds'))
evap <- evap[order(grid_id, date)]
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

cat("Grid cells:", nrow(evap_grid), "\n")
cat("Evap rows:", nrow(evap), "\n")

#===============================================================================
# 2. PENTAD STANDARDISATION
#===============================================================================
cat("\nCreating pentad statistics...\n")
pentads <- copy(evap)
pentads[, pentad := ceiling((yday(date) - leap_year(year(date)) * (yday(date) > 59)) / 5)]
pentads[, std_value := (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
        by = .(pentad, grid_id)]
pentads[, pentad_std_q95 := quantile(std_value, EXTREMES_THRES, na.rm = TRUE), by = grid_id]
pentads[, pentad_std_q80 := quantile(std_value, LOW_THRES, na.rm = TRUE), by = grid_id]
gc()

# Quantile-regression thresholds (non-stationarity correction)
cat("Computing quantile-regression thresholds (this may take a while)...\n")
pentads[, pentad_median_qr := tryCatch(
  rq(std_value ~ date, tau = 0.5)$fitted,
  error = function(e) rep(median(std_value), .N)
), by = .(pentad, grid_id)]

pentads[, pentad_std_q95_qr := tryCatch(
  rq(std_value ~ date, tau = EXTREMES_THRES)$fitted,
  error = function(e) rep(quantile(std_value, EXTREMES_THRES), .N)
), by = .(pentad, grid_id)]
pentads[, value := NULL]
gc()

cat("Pentad statistics complete.\n")

#===============================================================================
# 3. EVENT IDENTIFICATION - DEFINITION 1: Mean / Q95
#===============================================================================
cat("\nIdentifying events: Mean/Q95 definition...\n")
exeves <- merge(evap, pentads[, .(grid_id, date, std_value, pentad_median_qr,
                                   pentad_std_q80, pentad_std_q95, pentad_std_q95_qr)],
                all.x = TRUE, by = c("grid_id", "date"))
gc()

exeves[, evap_event := FALSE]
exeves[, value_above_low_thres := FALSE]
exeves[, extreme := FALSE]
exeves[std_value > 0, value_above_low_thres := TRUE]
exeves[std_value > pentad_std_q95, extreme := TRUE]
exeves[, above_low_thres_id := rleid(value_above_low_thres), by = grid_id]
exeves[, extreme_id := rleid(extreme), by = grid_id]

exeves[extreme == TRUE, evap_event := TRUE, .(grid_id, above_low_thres_id)]
above_low_thres_ids_with_extreme <- exeves[extreme == TRUE, above_low_thres_id]
exeves[above_low_thres_id %in% above_low_thres_ids_with_extreme, evap_event := TRUE]
exeves[, event_id := rleid(evap_event), .(grid_id)]
exeves[evap_event != TRUE, event_id := NA]
exeves[extreme != TRUE, extreme_id := NA]

# Add period and season
exeves[, period := ordered('up_to_2001')]
exeves[date > END_PERIOD_1, period := ordered('after_2001')]

exeves[month(date) < 4,  season := ordered("JFM")]
exeves[month(date) >= 4  & month(date) < 7,  season := ordered("AMJ")]
exeves[month(date) >= 7  & month(date) < 10, season := ordered("JAS")]
exeves[month(date) >= 10, season := ordered("OND")]
gc()

#===============================================================================
# 4. EVENT IDENTIFICATION - DEFINITION 2: Median (QR) / Q95 (QR)
#===============================================================================
cat("Identifying events: Median/Q95* definition (quantile regression)...\n")
exeves_qr <- merge(evap, pentads[, .(grid_id, date, std_value, pentad_median_qr,
                                      pentad_std_q95_qr)],
                   all.x = TRUE, by = c("grid_id", "date"))
exeves_qr[, evap_event := FALSE]
exeves_qr[, value_above_low_thres := FALSE]
exeves_qr[, extreme := FALSE]
exeves_qr[std_value > pentad_median_qr, value_above_low_thres := TRUE]
exeves_qr[std_value > pentad_std_q95_qr, extreme := TRUE]
exeves_qr[, above_low_thres_id := rleid(value_above_low_thres)]
exeves_qr[, extreme_qr_id := rleid(extreme), .(grid_id)]

exeves_qr[extreme == TRUE, evap_event := TRUE, .(grid_id, above_low_thres_id)]
above_low_thres_ids_with_extreme <- exeves_qr[extreme == TRUE, above_low_thres_id]
exeves_qr[above_low_thres_id %in% above_low_thres_ids_with_extreme, evap_event := TRUE]
exeves_qr[, event_qr_id := rleid(evap_event), .(grid_id)]
exeves_qr[evap_event != TRUE, event_qr_id := NA]
exeves_qr[extreme != TRUE, extreme_qr_id := NA]
gc()

#===============================================================================
# 5. EVENT IDENTIFICATION - DEFINITION 3: Q80 / Q95
#===============================================================================
cat("Identifying events: Q80/Q95 definition...\n")
exeves_80_95 <- merge(evap, pentads[, .(grid_id, date, std_value, pentad_median_qr,
                                         pentad_std_q80, pentad_std_q95, pentad_std_q95_qr)],
                      all.x = TRUE, by = c("grid_id", "date"))

exeves_80_95[, evap_event := FALSE]
exeves_80_95[, value_above_low_thres := FALSE]
exeves_80_95[, extreme := FALSE]
exeves_80_95[std_value > pentad_std_q80, value_above_low_thres := TRUE]
exeves_80_95[std_value > pentad_std_q95, extreme := TRUE]
exeves_80_95[, above_low_thres_id := rleid(value_above_low_thres)]
exeves_80_95[, extreme_id := rleid(extreme), .(grid_id)]

exeves_80_95[extreme == TRUE, evap_event := TRUE, .(grid_id, above_low_thres_id)]
above_low_thres_ids_with_extreme <- exeves_80_95[extreme == TRUE, above_low_thres_id]
exeves_80_95[above_low_thres_id %in% above_low_thres_ids_with_extreme, evap_event := TRUE]
exeves_80_95[, event_80_95_id := rleid(evap_event), .(grid_id)]
exeves_80_95[evap_event != TRUE, event_80_95_id := NA]
exeves_80_95[extreme != TRUE, extreme_id := NA]
gc()

#===============================================================================
# 6. EVENT IDENTIFICATION - DEFINITION 4: Q80 only
#===============================================================================
cat("Identifying events: Q80 definition...\n")
exeves_80 <- merge(evap, pentads[, .(grid_id, date, std_value, pentad_median_qr,
                                      pentad_std_q80)],
                   all.x = TRUE, by = c("grid_id", "date"))
exeves_80[, evap_event := FALSE]
exeves_80[std_value > pentad_std_q80, evap_event := TRUE]
exeves_80[, event_80_id := rleid(evap_event), .(grid_id)]
exeves_80[evap_event != TRUE, event_80_id := NA]
gc()

#===============================================================================
# 7. MERGE ALL DEFINITIONS
#===============================================================================
cat("Merging all definitions...\n")
exeves <- merge(exeves[, .(grid_id, date, season, period, value, std_value,
                            event_id, extreme_id)],
                exeves_qr[, .(grid_id, date, event_qr_id, extreme_qr_id)],
                by = c('grid_id', 'date'))
exeves <- merge(exeves,
                exeves_80_95[, .(grid_id, date, event_80_95_id)],
                by = c('grid_id', 'date'))
exeves <- merge(exeves,
                exeves_80[, .(grid_id, date, event_80_id)],
                by = c('grid_id', 'date'))

rm(exeves_qr, exeves_80_95, exeves_80); gc()

#===============================================================================
# 8. SAVE
#===============================================================================
cat("Saving results...\n")
saveRDS(pentads, paste0(PATH_OUTPUT_DATA, 'pentads_std_', region, '.rds'))
saveRDS(evap_grid, paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))
saveRDS(exeves, paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
cat("  Saved:", paste0('exeves_std_', region, '.rds'), "\n")

#===============================================================================
# 9. SUMMARY
#===============================================================================
cat("\n=== PREPROCESSING SUMMARY ===\n")
cat("Grid cells     :", exeves[, uniqueN(grid_id)], "\n")
cat("Date range     :", as.character(min(exeves$date)), "to", as.character(max(exeves$date)), "\n")
cat("Events/cell (Q80/Q95) :", round(mean(exeves[!is.na(event_80_95_id),
  uniqueN(event_80_95_id), grid_id]$V1), 1), "\n")
cat("Events/cell (Mean/Q95):", round(mean(exeves[!is.na(event_id),
  uniqueN(event_id), grid_id]$V1), 1), "\n")
cat("Events/cell (QR)      :", round(mean(exeves[!is.na(event_qr_id),
  uniqueN(event_qr_id), grid_id]$V1), 1), "\n")
cat("Events/cell (Q80)     :", round(mean(exeves[!is.na(event_80_id),
  uniqueN(event_80_id), grid_id]$V1), 1), "\n")

rm(evap, pentads, exeves); gc()
cat("\nPreprocessing complete!\n")
