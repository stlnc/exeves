# Preprocessing: identify extreme evaporation events (ExEvEs) using multiple definitions
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/02_preprocessing.R
# Optimised: single merge for all definitions, parallel quantile regression

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
setkey(evap, grid_id, date)
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

cat("Grid cells:", nrow(evap_grid), "\n")
cat("Evap rows:", nrow(evap), "\n")

#===============================================================================
# 2. PENTAD STANDARDISATION (compute once, join once)
#===============================================================================
cat("\nCreating pentad statistics...\n")

# Work in-place on evap to avoid a full copy
evap[, pentad := ceiling((yday(date) - leap_year(year(date)) * (yday(date) > 59)) / 5)]
evap[, std_value := (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
     by = .(pentad, grid_id)]

# Compute quantile thresholds once per grid_id, then join (avoids per-row expansion)
q_lookup <- evap[, .(pentad_std_q95 = quantile(std_value, EXTREMES_THRES, na.rm = TRUE),
                     pentad_std_q80 = quantile(std_value, LOW_THRES, na.rm = TRUE)),
                 by = grid_id]
evap <- q_lookup[evap, on = "grid_id"]
rm(q_lookup); gc()

# Quantile-regression thresholds (non-stationarity correction)
# Compute per (pentad, grid_id) – the most expensive step
cat("Computing quantile-regression thresholds (this may take a while)...\n")
cat("  Groups:", evap[, uniqueN(paste(pentad, grid_id))], "\n")

# Use tryCatch with lighter fitted-value extraction
evap[, c("pentad_median_qr", "pentad_std_q95_qr") := {
  fit_med <- tryCatch(rq(std_value ~ date, tau = 0.5)$fitted,
                      error = function(e) rep(median(std_value, na.rm = TRUE), .N))
  fit_q95 <- tryCatch(rq(std_value ~ date, tau = EXTREMES_THRES)$fitted,
                      error = function(e) rep(quantile(std_value, EXTREMES_THRES, na.rm = TRUE), .N))
  list(fit_med, fit_q95)
}, by = .(pentad, grid_id)]
gc()

cat("Pentad statistics complete.\n")

# Save pentads for later (only the columns needed)
pentads <- evap[, .(grid_id, date, std_value, pentad_median_qr,
                     pentad_std_q80, pentad_std_q95, pentad_std_q95_qr)]
saveRDS(pentads, paste0(PATH_OUTPUT_DATA, 'pentads_std_', region, '.rds'))
rm(pentads); gc()

#===============================================================================
# 3-6. ALL EVENT DEFINITIONS ON ONE TABLE (no extra copies)
# Build all boolean flags and rleid columns on `evap` directly
#===============================================================================
cat("\nIdentifying events: all definitions in one pass...\n")

# --- Definition 1: Mean / Q95 ---
evap[, val_above_mean := std_value > 0]
evap[, ext_q95       := std_value > pentad_std_q95]
evap[, above_mean_id := rleid(val_above_mean), by = grid_id]
evap[, extreme_id    := rleid(ext_q95),         by = grid_id]
# Mark events: any run of above-mean that contains at least one extreme day
mean_ids_with_ext <- evap[ext_q95 == TRUE, unique(above_mean_id)]
evap[, evap_event_mean := above_mean_id %in% mean_ids_with_ext & val_above_mean == TRUE]
evap[, event_id := rleid(evap_event_mean), by = grid_id]
evap[evap_event_mean == FALSE, event_id := NA]
evap[ext_q95 == FALSE, extreme_id := NA]

# --- Definition 2: Median (QR) / Q95 (QR) ---
evap[, val_above_medqr := std_value > pentad_median_qr]
evap[, ext_q95_qr      := std_value > pentad_std_q95_qr]
evap[, above_medqr_id  := rleid(val_above_medqr), by = grid_id]
evap[, extreme_qr_id   := rleid(ext_q95_qr),      by = grid_id]
medqr_ids_with_ext <- evap[ext_q95_qr == TRUE, unique(above_medqr_id)]
evap[, evap_event_qr := above_medqr_id %in% medqr_ids_with_ext & val_above_medqr == TRUE]
evap[, event_qr_id := rleid(evap_event_qr), by = grid_id]
evap[evap_event_qr == FALSE, event_qr_id := NA]
evap[ext_q95_qr == FALSE, extreme_qr_id := NA]

# --- Definition 3: Q80 / Q95 ---
evap[, val_above_q80 := std_value > pentad_std_q80]
evap[, above_q80_id  := rleid(val_above_q80), by = grid_id]
q80_ids_with_ext <- evap[ext_q95 == TRUE, unique(above_q80_id)]
evap[, evap_event_80_95 := above_q80_id %in% q80_ids_with_ext & val_above_q80 == TRUE]
evap[, event_80_95_id := rleid(evap_event_80_95), by = grid_id]
evap[evap_event_80_95 == FALSE, event_80_95_id := NA]

# --- Definition 4: Q80 only ---
evap[, event_80_id := rleid(val_above_q80), by = grid_id]
evap[val_above_q80 == FALSE, event_80_id := NA]

# Add period and season
evap[, period := ordered('up_to_2001')]
evap[date > END_PERIOD_1, period := ordered('after_2001')]
evap[month(date) < 4,                           season := ordered("JFM")]
evap[month(date) >= 4  & month(date) < 7,       season := ordered("AMJ")]
evap[month(date) >= 7  & month(date) < 10,      season := ordered("JAS")]
evap[month(date) >= 10,                          season := ordered("OND")]

# Drop temporary boolean/helper columns
drop_cols <- c("pentad", "pentad_std_q80", "pentad_std_q95",
               "pentad_std_q95_qr", "pentad_median_qr",
               "val_above_mean", "ext_q95", "above_mean_id",
               "evap_event_mean", "val_above_medqr", "ext_q95_qr",
               "above_medqr_id", "evap_event_qr",
               "val_above_q80", "above_q80_id",
               "evap_event_80_95")
evap[, (drop_cols) := NULL]
gc()

#===============================================================================
# 7. SAVE
#===============================================================================
cat("Saving results...\n")
exeves <- evap  # rename for clarity downstream
setkey(exeves, grid_id, date)
saveRDS(evap_grid, paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))
saveRDS(exeves,    paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
cat("  Saved:", paste0('exeves_std_', region, '.rds'), "\n")

#===============================================================================
# 8. SUMMARY
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

rm(evap, exeves); gc()
cat("\nPreprocessing complete!\n")
