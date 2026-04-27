library(data.table)
library(terra)
library(lubridate)

aet_file    <- "gleam_e_mm_med_198001_202412_025_daily.nc"
precip_file <- "mswep-v2-8_tp_mm_med_197901_202012_025_daily.nc"

# ── 1. Load rasters ─────────────────────────────────────────────────────────
message("Loading rasters...")
R_aet    <- rast(aet_file)
R_precip <- rast(precip_file)

# ── 2. Verify spatial grids are identical ───────────────────────────────────
if (!compareGeom(R_aet, R_precip, stopOnError = FALSE, lyrs = FALSE,
                 crs = FALSE, ext = TRUE, rowcol = TRUE)) {
  stop("Spatial grids differ between AET and precipitation files.")
}
message(sprintf("Grids match: %d cells, %.4f deg resolution.",
                ncell(R_aet), res(R_aet)[1]))

# ── 3. Align to common time period ──────────────────────────────────────────
dates_aet    <- as.IDate(time(R_aet))
dates_precip <- as.IDate(time(R_precip))
common_dates <- sort(as.IDate(intersect(dates_aet, dates_precip)))

message(sprintf("Common period: %s to %s  (%d days)",
                min(common_dates), max(common_dates), length(common_dates)))

R_aet_sub    <- R_aet[[which(dates_aet    %in% common_dates)]]
R_precip_sub <- R_precip[[which(dates_precip %in% common_dates)]]

stopifnot(nlyr(R_aet_sub) == nlyr(R_precip_sub),
          nlyr(R_aet_sub) == length(common_dates))
message(sprintf("Both rasters aligned to %d layers.", nlyr(R_aet_sub)))

# ── 4. Detect ExEvEs (AET) ──────────────────────────────────────────────────
message("Building AET data.table...")
nc <- ncell(R_aet_sub)
nl <- nlyr(R_aet_sub)

DT <- data.table(
  cell  = rep(seq_len(nc), times = nl),
  date  = rep(common_dates, each = nc),
  val_1 = as.vector(values(R_aet_sub))
)
DT <- DT[!is.na(val_1)]
setorder(DT, cell, date)

message("Computing pentads and z-scores...")
DT[, pentad := ceiling(
  (yday(date) - as.integer(leap_year(date)) * (yday(date) > 59L)) / 5
)]
DT[pentad > 73L, pentad := 73L]

DT[, `:=`(mu_1 = mean(val_1, na.rm = TRUE),
           sd_1 = sd(val_1,   na.rm = TRUE)), by = .(cell, pentad)]
DT[, z_1 := (val_1 - mu_1) / sd_1]

DT[, `:=`(z1_thresh  = quantile(z_1, 0.80, na.rm = TRUE),
           z1_extreme = quantile(z_1, 0.95, na.rm = TRUE)), by = cell]

DT[, above_1   := z_1 >= z1_thresh]
DT[, extreme_1 := z_1 >= z1_extreme]

DT[, run_1 := rleid(above_1), by = cell]
DT[, `:=`(run_len_1 = .N,
           has_ext_1 = any(extreme_1)), by = .(cell, run_1)]
DT[, high_event_1 := above_1 & run_len_1 >= 2 & has_ext_1]
DT[, high_event_1_id := ifelse(high_event_1, rleid(high_event_1), NA_integer_),
   by = cell]
DT[, c("run_1", "run_len_1", "has_ext_1") := NULL]

n_exeve_days  <- sum(DT$high_event_1)
n_exeve_cells <- DT[high_event_1 == TRUE, uniqueN(cell)]
message(sprintf("ExEvE detection done: %d event-days across %d grid cells.",
                n_exeve_days, n_exeve_cells))

# ── 5. Build precipitation data.table ───────────────────────────────────────
message("Building precipitation data.table...")
DT_precip <- data.table(
  cell   = rep(seq_len(nc), times = nl),
  date   = rep(common_dates, each = nc),
  precip = as.vector(values(R_precip_sub))
)
DT_precip <- DT_precip[!is.na(precip)]
setorder(DT_precip, cell, date)

# ── 6. Compute q95 per grid cell ────────────────────────────────────────────
message("Computing q95 threshold per grid cell...")
DT_precip[, q95 := quantile(precip, 0.95, na.rm = TRUE), by = cell]
DT_precip[, is_peak := precip >= q95]

# ── 7. Flag precipitation days that fall inside an ExEvE ────────────────────
message("Matching peak precipitation days with ExEvE days...")
exeve_days <- DT[high_event_1 == TRUE, .(cell, date)]
setkey(exeve_days, cell, date)
setkey(DT_precip,  cell, date)

DT_precip[, is_exeve := FALSE]
DT_precip[exeve_days, on = .(cell, date), is_exeve := TRUE]

# ── 8. Aggregate per grid cell ──────────────────────────────────────────────
message("Aggregating results...")
result <- DT_precip[, .(
  n_total_days      = .N,
  n_exeve_days      = sum(is_exeve),
  n_peak_days       = sum(is_peak),
  n_peak_in_exeve   = sum(is_peak & is_exeve),
  q95_threshold_mm  = first(q95)
), by = cell]

result[, pct_peak_in_exeve := round(100 * n_peak_in_exeve / n_peak_days, 2)]

# Add coordinates
coords <- as.data.table(terra::xyFromCell(R_precip_sub, result$cell))
setnames(coords, c("lon", "lat"))
coords[, cell := result$cell]
result <- merge(result, coords, by = "cell", all.x = TRUE)
setcolorder(result, c("cell", "lon", "lat",
                      setdiff(names(result), c("cell", "lon", "lat"))))

# ── 9. Save ──────────────────────────────────────────────────────────────────
fwrite(result, "precip_peaks_during_ExEvEs.csv")
message("Saved to precip_peaks_during_ExEvEs.csv")

message("\n── Summary ────────────────────────────────────────────")
message(sprintf("  Grid cells analysed          : %d", nrow(result)))
message(sprintf("  Cells with any ExEvE days    : %d", sum(result$n_exeve_days > 0)))
message(sprintf("  Total peak precip days (q95) : %d", sum(result$n_peak_days)))
message(sprintf("  Peak precip days in ExEvEs   : %d", sum(result$n_peak_in_exeve)))
message(sprintf("  Mean %% peak days in ExEvEs   : %.2f%%",
                mean(result$pct_peak_in_exeve[result$n_exeve_days > 0], na.rm = TRUE)))
