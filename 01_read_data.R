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
# HELPER: Read NetCDF → data.table  (vectorised – no pixel loop)
#===============================================================================
read_nc_to_dt <- function(nc_path, var_name = NULL) {
  cat("  Opening:", nc_path, "\n")
  nc <- nc_open(nc_path)
  
  # Auto-detect variable name if not given
  if (is.null(var_name)) {
    var_names <- names(nc$var)
    var_name <- var_names[!var_names %in% c("time_bnds", "lat_bnds", "lon_bnds",
                                             "crs", "spatial_ref")][1]
  }
  cat("  Variable:", var_name, "\n")
  
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
  
  n_lon  <- length(lon)
  n_lat  <- length(lat)
  n_time <- length(time_raw)
  cat("  Dimensions: lon =", n_lon, ", lat =", n_lat, ", time =", n_time, "\n")
  cat("  Date range:", as.character(min(dates)), "to", as.character(max(dates)), "\n")
  
  # --- Subset time dimension BEFORE reading (if file supports it) ---
  # Find which time steps fall within our analysis period
  t_idx <- which(dates >= START_PERIOD_1 & dates <= END_PERIOD_2)
  if (length(t_idx) > 0 && length(t_idx) < n_time) {
    cat("  Subsetting time in-file:", length(t_idx), "of", n_time, "steps\n")
    vals  <- ncvar_get(nc, var_name,
                       start = c(1, 1, min(t_idx)),
                       count = c(-1, -1, length(t_idx)))
    dates <- dates[t_idx]
    n_time <- length(t_idx)
  } else {
    vals <- ncvar_get(nc, var_name)
  }
  nc_close(nc)
  
  # --- Vectorised reshape: 3-D array → data.table in one shot ---
  # Identify valid (non-all-NA) pixels via fast colSums on reshaped matrix
  dim(vals) <- c(n_lon * n_lat, n_time)           # pixels × time
  pixel_has_data <- rowSums(!is.na(vals)) > 0L
  valid_idx <- which(pixel_has_data)
  n_valid   <- length(valid_idx)
  cat("  Valid pixels:", n_valid, "\n")
  
  # Build coordinate vectors for valid pixels
  lon_idx <- ((valid_idx - 1L) %% n_lon) + 1L
  lat_idx <- ((valid_idx - 1L) %/% n_lon) + 1L
  
  # Extract valid sub-matrix and build data.table in one allocation
  vals_valid <- vals[valid_idx, , drop = FALSE]     # n_valid × n_time
  rm(vals); gc()
  
  dt <- data.table(
    lon   = rep(lon[lon_idx], each = n_time),
    lat   = rep(lat[lat_idx], each = n_time),
    date  = rep(dates, times = n_valid),
    value = as.vector(t(vals_valid))                # row-major unroll
  )
  rm(vals_valid); gc()
  
  # Drop any remaining individual NAs
  dt <- dt[!is.na(value)]
  
  return(dt)
}

#===============================================================================
# 1. READ EVAPORATION (GLEAM) — time already subsetted in-file
#===============================================================================
cat("\n=== Reading GLEAM evaporation ===\n")
evap_raw <- read_nc_to_dt(EVAP_NC_FILE)
cat("  Rows:", nrow(evap_raw), "\n")

#===============================================================================
# 2. READ PRECIPITATION (MSWEP) — time already subsetted in-file
#===============================================================================
cat("\n=== Reading MSWEP precipitation ===\n")
prec_raw <- read_nc_to_dt(PREC_NC_FILE)
cat("  Rows:", nrow(prec_raw), "\n")

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

# Keep only common grid points – use keyed semi-join
common_coords <- fintersect(evap_coords, prec_coords)
cat("  Common grid points:", nrow(common_coords), "\n")

setkey(evap_raw, lon, lat)
setkey(prec_raw, lon, lat)
evap_raw <- evap_raw[common_coords, nomatch = 0]
prec_raw <- prec_raw[common_coords, nomatch = 0]

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
evap_raw[, grid_id := .GRP, by = .(lat, lon)]
evap_grid <- unique(evap_raw[, .(lon, lat, grid_id)])
GRID_CELL_N <- nrow(evap_grid)
cat("  Number of grid cells:", GRID_CELL_N, "\n")

# Save grid
saveRDS(evap_grid, paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

# Evaporation: keep only needed columns, set key, save
evap <- evap_raw[, .(grid_id, date, value)]
setkey(evap, grid_id, date)
saveRDS(evap, paste0(PATH_OUTPUT_DATA, region, '_evap_grid.rds'))
cat("  Saved:", paste0(region, '_evap_grid.rds'), "\n")

# Precipitation: join grid IDs via keyed lookup (no extra copy)
setkey(evap_grid, lon, lat)
prec_raw[, `:=`(lon = round(lon, 4), lat = round(lat, 4))]
setkey(prec_raw, lon, lat)
prec <- evap_grid[prec_raw, nomatch = 0][, .(grid_id, date, value)]
setkey(prec, grid_id, date)

# Add extreme precipitation flag (90th percentile per grid cell)
q90_lookup <- prec[, .(q90 = quantile(value, 0.9, na.rm = TRUE)), grid_id]
prec <- q90_lookup[prec, on = "grid_id"]
prec[value > q90, extreme_prec := TRUE]
prec[, q90 := NULL]

saveRDS(prec, paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))
cat("  Saved:", paste0(region, '_prec_grid.rds'), "\n")

# Save GRID_CELL_N for later scripts
save(GRID_CELL_N, file = paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

# Clean up
rm(evap_raw, prec_raw, prec, common_coords, evap_coords, prec_coords, q90_lookup); gc()

cat("\n=== Data reading complete ===\n")
cat("  Grid cells:", GRID_CELL_N, "\n")
cat("  Common period:", as.character(common_dates_start), "to",
    as.character(common_dates_end), "\n")
