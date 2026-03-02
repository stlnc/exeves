# Read and align evaporation (GLEAM) and precipitation (MSWEP) NetCDF data
# Both datasets share 0.25° resolution; subset to common analysis period
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/01_crop_data.R

library(data.table)
library(ncdf4)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

#===============================================================================
# HELPER: Read NetCDF → data.table
#===============================================================================
read_nc_to_dt <- function(nc_path, var_name = NULL) {
  cat("  Opening:", nc_path, "\n")
  nc <- nc_open(nc_path)
  
  # Auto-detect variable name if not given
  if (is.null(var_name)) {
    var_names <- names(nc$var)
    # Skip coordinate/bounds variables
    var_name <- var_names[!var_names %in% c("time_bnds", "lat_bnds", "lon_bnds",
                                             "crs", "spatial_ref")][1]
  }
  cat("  Variable:", var_name, "\n")
  
  vals <- ncvar_get(nc, var_name)
  lon  <- ncvar_get(nc, "lon")
  lat  <- ncvar_get(nc, "lat")
  time_raw <- ncvar_get(nc, "time")
  
  # Convert time to dates
  time_units <- ncatt_get(nc, "time", "units")$value
  time_origin <- as.Date(sub(".*since ", "", time_units))
  if (grepl("hours", time_units)) {
    dates <- time_origin + time_raw / 24
  } else {
    dates <- time_origin + time_raw
  }
  nc_close(nc)
  
  cat("  Dimensions: lon =", length(lon), ", lat =", length(lat),
      ", time =", length(dates), "\n")
  cat("  Date range:", as.character(min(dates)), "to", as.character(max(dates)), "\n")
  
  # Build data.table
  dt_list <- vector("list", length(lon) * length(lat))
  k <- 0L
  for (i in seq_along(lon)) {
    for (j in seq_along(lat)) {
      ts <- vals[i, j, ]
      if (!all(is.na(ts))) {
        k <- k + 1L
        dt_list[[k]] <- data.table(lon = lon[i], lat = lat[j],
                                   date = dates, value = ts)
      }
    }
  }
  rm(vals); gc()
  dt <- rbindlist(dt_list[seq_len(k)])
  rm(dt_list); gc()
  dt <- dt[!is.na(value)]
  dt <- dt[order(lon, lat, date)]
  cat("  Valid pixels:", dt[, uniqueN(paste(lon, lat))], "\n")
  return(dt)
}

#===============================================================================
# 1. READ EVAPORATION (GLEAM)
#===============================================================================
cat("\n=== Reading GLEAM evaporation ===\n")
evap_raw <- read_nc_to_dt(EVAP_NC_FILE)

# Subset to analysis period
evap_raw <- evap_raw[date >= START_PERIOD_1 & date <= END_PERIOD_2]
cat("  After period filter:", nrow(evap_raw), "rows\n")
cat("  Date range:", as.character(min(evap_raw$date)), "to",
    as.character(max(evap_raw$date)), "\n")

#===============================================================================
# 2. READ PRECIPITATION (MSWEP)
#===============================================================================
cat("\n=== Reading MSWEP precipitation ===\n")
prec_raw <- read_nc_to_dt(PREC_NC_FILE)

# Subset to analysis period
prec_raw <- prec_raw[date >= START_PERIOD_1 & date <= END_PERIOD_2]
cat("  After period filter:", nrow(prec_raw), "rows\n")
cat("  Date range:", as.character(min(prec_raw$date)), "to",
    as.character(max(prec_raw$date)), "\n")

#===============================================================================
# 3. ALIGN SPATIAL GRIDS
# Both are 0.25° but ensure exact lon/lat matching
#===============================================================================
cat("\n=== Aligning grids ===\n")

# Round coordinates to avoid floating-point mismatches
evap_raw[, `:=`(lon = round(lon, 4), lat = round(lat, 4))]
prec_raw[, `:=`(lon = round(lon, 4), lat = round(lat, 4))]

evap_coords <- unique(evap_raw[, .(lon, lat)])
prec_coords <- unique(prec_raw[, .(lon, lat)])
cat("  Evap grid points:", nrow(evap_coords), "\n")
cat("  Prec grid points:", nrow(prec_coords), "\n")

# Keep only common grid points
common_coords <- merge(evap_coords, prec_coords, by = c("lon", "lat"))
cat("  Common grid points:", nrow(common_coords), "\n")

evap_raw <- evap_raw[common_coords, on = .(lon, lat), nomatch = 0]
prec_raw <- prec_raw[common_coords, on = .(lon, lat), nomatch = 0]

#===============================================================================
# 4. ENSURE IDENTICAL DATE RANGE
#===============================================================================
cat("\n=== Aligning dates ===\n")
common_dates_start <- max(min(evap_raw$date), min(prec_raw$date))
common_dates_end   <- min(max(evap_raw$date), max(prec_raw$date))
cat("  Common date range:", as.character(common_dates_start), "to",
    as.character(common_dates_end), "\n")

evap_raw <- evap_raw[date >= common_dates_start & date <= common_dates_end]
prec_raw <- prec_raw[date >= common_dates_start & date <= common_dates_end]

cat("  Evap rows:", nrow(evap_raw), "\n")
cat("  Prec rows:", nrow(prec_raw), "\n")

#===============================================================================
# 5. ASSIGN GRID IDs AND SAVE
#===============================================================================
cat("\n=== Assigning grid IDs and saving ===\n")

# Grid IDs based on evaporation grid (same as prec since aligned)
evap_raw <- evap_raw[order(lon, lat)]
evap_raw[, grid_id := .GRP, by = list(lat, lon)]
evap_grid <- unique(evap_raw[, .(lon, lat, grid_id)])
GRID_CELL_N <- nrow(evap_grid)
cat("  Number of grid cells:", GRID_CELL_N, "\n")

# Save grid
saveRDS(evap_grid, paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

# Evaporation: save with grid_id
evap <- evap_raw[, .(grid_id, date, value)]
saveRDS(evap, paste0(PATH_OUTPUT_DATA, region, '_evap_grid.rds'))
cat("  Saved:", paste0(region, '_evap_grid.rds'), "\n")

# Precipitation: merge grid IDs and save
prec_raw[, `:=`(lon = round(lon, 4), lat = round(lat, 4))]
prec <- merge(prec_raw, evap_grid, by = c("lon", "lat"))
prec <- prec[, .(grid_id, date, value)]

# Add extreme precipitation flag (90th percentile per grid cell)
prec[, q90 := quantile(value, 0.9, na.rm = TRUE), grid_id]
prec[value > q90, extreme_prec := TRUE][, q90 := NULL]

saveRDS(prec, paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))
cat("  Saved:", paste0(region, '_prec_grid.rds'), "\n")

# Save GRID_CELL_N for later scripts
save(GRID_CELL_N, file = paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

# Clean up
rm(evap_raw, prec_raw, prec, common_coords, evap_coords, prec_coords); gc()

cat("\n=== Data reading complete ===\n")
cat("  Grid cells:", GRID_CELL_N, "\n")
cat("  Common period:", as.character(common_dates_start), "to",
    as.character(common_dates_end), "\n")
