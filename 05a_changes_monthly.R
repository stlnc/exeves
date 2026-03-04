# Monthly change analysis: evaporation and precipitation by ExEvE/non-ExEvE conditions
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/05a_changes_monthly.R
# Note: radiation data not available; using evaporation + precipitation only

library(data.table)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec   <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))

# Merge with precipitation – keyed join
setkey(exeves, grid_id, date)
setkey(prec, grid_id, date)
exeves_changes <- prec[, .(grid_id, date, prec = value)][exeves, on = .(grid_id, date), nomatch = 0]
exeves_changes <- exeves_changes[, .(grid_id, period, month = month(date, label = TRUE),
                                      event_80_95_id, evap = value, prec)]

exeves_changes[, conditions := ordered('ExEvE')]
exeves_changes[is.na(event_80_95_id), conditions := ordered('non-ExEvE')]
exeves_changes[, event_80_95_id := NULL]

# Monthly sums
exeves_changes_summary <- exeves_changes[, .(evap = sum(evap),
                                              prec = sum(prec)),
                                          by = .(grid_id, month, period, conditions)]

# Monthly means
exeves_changes_mean_summary <- exeves_changes[, .(evap = mean(evap),
                                                   prec = mean(prec)),
                                               by = .(grid_id, month, period, conditions)]

saveRDS(exeves_changes_summary,      file = paste0(PATH_OUTPUT, 'monthly_changes.rds'))
saveRDS(exeves_changes_mean_summary, file = paste0(PATH_OUTPUT, 'monthly_changes_mean.rds'))

cat("Monthly change data saved.\n")
rm(exeves, prec, exeves_changes, exeves_changes_summary, exeves_changes_mean_summary); gc()
