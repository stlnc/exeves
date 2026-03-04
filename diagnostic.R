# diagnostic.R — Run on HPC to gather debugging info for 05b and 07c
# Usage:  Rscript diagnostic.R  (or source from R session)
# Output: diagnostics.txt in the project root

library(data.table)
library(lubridate)

load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

sink("diagnostics.txt", split = TRUE)

cat("============================================================\n")
cat("  ExEvEs DIAGNOSTIC REPORT\n")
cat("  Generated:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ── 1. ENVIRONMENT ──────────────────────────────────────────────
cat("=== 1. ENVIRONMENT ===\n")
cat("Region          :", region, "\n")
cat("GRID_CELL_N     :", GRID_CELL_N, "\n")
cat("PERIOD_LENGTH   :", PERIOD_LENGTH, "\n")
cat("END_PERIOD_1    :", as.character(END_PERIOD_1), "\n")
cat("END_PERIOD_2    :", as.character(END_PERIOD_2), "\n")
cat("EXTREMES_THRES  :", EXTREMES_THRES, "\n")
cat("LOW_THRES       :", LOW_THRES, "\n\n")

# ── 2. EXEVES TABLE ────────────────────────────────────────────
cat("=== 2. EXEVES TABLE ===\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
cat("Dimensions   :", nrow(exeves), "rows x", ncol(exeves), "cols\n")
cat("Columns      :", paste(colnames(exeves), collapse = ", "), "\n")
cat("Grid cells   :", exeves[, uniqueN(grid_id)], "\n")
cat("Date range   :", as.character(min(exeves$date)), "to",
    as.character(max(exeves$date)), "\n")
cat("Period counts:\n")
print(exeves[, .N, by = period])
cat("\nSample rows (first grid cell, first 5 dates):\n")
first_gid <- exeves[, min(grid_id)]
print(exeves[grid_id == first_gid][1:5])
cat("\n")

# ── 3. PREC TABLE ──────────────────────────────────────────────
cat("=== 3. PRECIPITATION TABLE ===\n")
prec <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))
cat("Dimensions   :", nrow(prec), "rows x", ncol(prec), "cols\n")
cat("Columns      :", paste(colnames(prec), collapse = ", "), "\n")
cat("Grid cells   :", prec[, uniqueN(grid_id)], "\n")
cat("Date range   :", as.character(min(prec$date)), "to",
    as.character(max(prec$date)), "\n")
cat("NA in value  :", sum(is.na(prec$value)), "\n")
cat("Zeros        :", sum(prec$value == 0, na.rm = TRUE), "\n")
cat("value summary:\n")
print(summary(prec$value))
cat("\n")

# ── 4. GRID ALIGNMENT ──────────────────────────────────────────
cat("=== 4. GRID ALIGNMENT ===\n")
eg <- unique(exeves$grid_id)
pg <- unique(prec$grid_id)
cat("Exeves grids :", length(eg), "\n")
cat("Prec grids   :", length(pg), "\n")
cat("Common grids :", length(intersect(eg, pg)), "\n")
cat("Exeves-only  :", length(setdiff(eg, pg)), "\n")
cat("Prec-only    :", length(setdiff(pg, eg)), "\n")

# Check date coverage per grid
edpg <- exeves[, .N, by = grid_id]
pdpg <- prec[, .N, by = grid_id]
cat("Exeves dates/grid : min=", min(edpg$N), " max=", max(edpg$N),
    " median=", median(edpg$N), "\n")
cat("Prec dates/grid   : min=", min(pdpg$N), " max=", max(pdpg$N),
    " median=", median(pdpg$N), "\n\n")

# ── 5. EVENT DEFINITIONS ───────────────────────────────────────
cat("=== 5. EVENT DEFINITIONS ===\n")
for (def in c("event_id", "event_qr_id", "event_80_95_id", "event_80_id")) {
  if (!def %in% colnames(exeves)) {
    cat(def, ": NOT PRESENT\n\n")
    next
  }
  sub <- exeves[!is.na(get(def))]
  durations <- sub[, .N, by = .(grid_id, get(def))]
  dur_dist  <- durations[, .N, by = .(duration = N)][order(duration)]

  cat(def, ":\n")
  cat("  Event-day rows   :", nrow(sub), "\n")
  cat("  Unique events    :", nrow(durations), "\n")
  cat("  Mean duration    :", round(mean(durations$N), 2), "\n")
  cat("  Median duration  :", median(durations$N), "\n")
  cat("  Duration distribution (top 15):\n")
  print(dur_dist[order(-N)][1:min(15, nrow(dur_dist))])
  cat("\n")
}

# ── 6. EXTREME_ID ──────────────────────────────────────────────
cat("=== 6. EXTREME_ID ===\n")
if ("extreme_id" %in% colnames(exeves)) {
  cat("Non-NA rows    :", sum(!is.na(exeves$extreme_id)), "\n")
  cat("Unique extremes:", exeves[!is.na(extreme_id),
                                 uniqueN(paste(grid_id, extreme_id))], "\n\n")
}

# ── 7. STD_VALUE & THRESHOLDS (single grid cell) ──────────────
cat("=== 7. STD_VALUE CHECK (grid_id =", first_gid, ") ===\n")
if ("std_value" %in% colnames(exeves)) {
  sg <- exeves[grid_id == first_gid]
  cat("std_value summary:\n")
  print(summary(sg$std_value))
  cat("Fraction > 0 (above mean)  :", round(mean(sg$std_value > 0, na.rm = TRUE), 4), "\n")
  # Check if pentad thresholds were saved
  pent_file <- paste0(PATH_OUTPUT_DATA, 'pentads_std_', region, '.rds')
  if (file.exists(pent_file)) {
    pent <- readRDS(pent_file)
    pg1 <- pent[grid_id == first_gid]
    cat("pentad_std_q80 (grid 1)    :", round(unique(pg1$pentad_std_q80)[1], 4), "\n")
    cat("pentad_std_q95 (grid 1)    :", round(unique(pg1$pentad_std_q95)[1], 4), "\n")
    cat("Fraction > q80             :", round(mean(pg1$std_value > pg1$pentad_std_q80, na.rm = TRUE), 4), "\n")
    cat("Fraction > q95             :", round(mean(pg1$std_value > pg1$pentad_std_q95, na.rm = TRUE), 4), "\n")
    rm(pg1, pent)
  } else {
    cat("pentads file not found, skipping threshold check\n")
  }
} else {
  cat("std_value column not found in exeves\n")
}
cat("\n")

# ── 8. 07c SIMULATION ─────────────────────────────────────────
cat("=== 8. 07c SIMULATION (Panel A/B) ===\n")
setkey(exeves, grid_id, date)
setkey(prec, grid_id, date)
exeves[prec, prec_val := i.value, on = .(grid_id, date)]

exeves_prec <- exeves[!is.na(event_80_95_id)]
exeves_prec[, event_duration := .N, by = .(event_80_95_id, grid_id)]

dur_counts <- exeves_prec[, .(N_events = uniqueN(paste(grid_id, event_80_95_id))),
                           by = event_duration][order(event_duration)]
cat("Event duration distribution (event_80_95_id):\n")
print(dur_counts)

MODAL_DUR <- dur_counts[which.max(N_events), event_duration]
cat("\nMODAL_DURATION :", MODAL_DUR, "\n")
cat("Events at modal:", dur_counts[event_duration == MODAL_DUR, N_events], "\n")

# Also show top-5 durations
cat("\nTop 5 most common durations:\n")
print(dur_counts[order(-N_events)][1:min(5, nrow(dur_counts))])
cat("\n")

# ── 9. 07c PANEL C SIMULATION ─────────────────────────────────
cat("=== 9. 07c PANEL C (Wet/Dry ratio) ===\n")
exeves[, prec_day := factor(fifelse(prec_val >= 1, "wet", "dry"),
                             levels = c("wet", "dry"))]
cat("prec_val NA count   :", sum(is.na(exeves$prec_val)), "\n")
cat("prec_day NA count   :", sum(is.na(exeves$prec_day)), "\n")
cat("Wet days (all)      :", sum(exeves$prec_day == "wet", na.rm = TRUE), "\n")
cat("Dry days (all)      :", sum(exeves$prec_day == "dry", na.rm = TRUE), "\n")
cat("Wet days (ExEvEs)   :",
    sum(exeves[!is.na(event_80_95_id)]$prec_day == "wet", na.rm = TRUE), "\n")
cat("Dry days (ExEvEs)   :",
    sum(exeves[!is.na(event_80_95_id)]$prec_day == "dry", na.rm = TRUE), "\n")

# Sample Panel C table for month=1
cat("\nPanel C data for month 1 (ExEvEs):\n")
exeves_prec_days <- exeves[!is.na(event_80_95_id) & !is.na(prec_val),
                            .N, .(period, prec_day, month = month(date))]
cat("  Counts:\n")
print(dcast(exeves_prec_days[month == 1], period ~ prec_day, value.var = "N"))
cat("  Ratios:\n")
totals <- exeves_prec_days[month == 1, .(total = sum(N)), by = period]
merged <- merge(exeves_prec_days[month == 1], totals, by = "period")
merged[, ratio := round(N / total, 4)]
print(merged[, .(period, prec_day, N, ratio)])
cat("\n")

# ── 10. 05b SIMULATION ────────────────────────────────────────
cat("=== 10. 05b SIMULATION (Spatial changes) ===\n")
evap_grid <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))
setkey(evap_grid, grid_id)
exeves[evap_grid, `:=`(lon = i.lon, lat = i.lat), on = "grid_id"]

# Check inner join effect
n_total <- nrow(exeves)
n_with_prec <- sum(!is.na(exeves$prec_val))
cat("Total rows       :", n_total, "\n")
cat("Rows with prec   :", n_with_prec, "\n")
cat("Rows without prec:", n_total - n_with_prec, "\n")
cat("Drop fraction    :", round(1 - n_with_prec / n_total, 4), "\n\n")

# Severity (event_id, ExEvEs) with inner-join filter
exeves_ij <- exeves[!is.na(prec_val)]
event_sev <- exeves_ij[!is.na(event_id), .(value = sum(value)), .(lon, lat, period)]
n_both <- event_sev[, .N, by = .(lon, lat)][N == 2, .N]
n_one  <- event_sev[, .N, by = .(lon, lat)][N == 1, .N]
cat("Severity (event_id):\n")
cat("  Pixels total     :", nrow(unique(event_sev[, .(lon, lat)])), "\n")
cat("  Pixels both periods:", n_both, "\n")
cat("  Pixels one period  :", n_one, "\n")

# Ratio distribution
event_sev2 <- event_sev[event_sev[, .N, by = .(lon, lat)][N == 2, .(lon, lat)],
                         on = .(lon, lat), nomatch = NULL]
event_sev2[, diff_value := diff(value), by = .(lon, lat)]
event_sev2[, ratio := 1 + diff_value / value, by = .(lon, lat, period)]
cat("  Ratio (up_to_2001 period):\n")
r <- event_sev2[period == "up_to_2001", ratio]
cat("    Min   :", round(min(r, na.rm = TRUE), 4), "\n")
cat("    Q25   :", round(quantile(r, 0.25, na.rm = TRUE), 4), "\n")
cat("    Median:", round(median(r, na.rm = TRUE), 4), "\n")
cat("    Q75   :", round(quantile(r, 0.75, na.rm = TRUE), 4), "\n")
cat("    Max   :", round(max(r, na.rm = TRUE), 4), "\n")
cat("    NaN   :", sum(is.nan(r)), "\n")
cat("    Inf   :", sum(is.infinite(r)), "\n\n")

# ── 11. SINGLE GRID CELL TRACE ────────────────────────────────
cat("=== 11. SINGLE GRID CELL TRACE (grid_id =", first_gid, ") ===\n")
g1 <- exeves[grid_id == first_gid]
cat("Total days     :", nrow(g1), "\n")
cat("event_id events:", g1[!is.na(event_id), uniqueN(event_id)], "\n")
cat("event_80_95_id :", g1[!is.na(event_80_95_id), uniqueN(event_80_95_id)], "\n")
cat("extreme_id     :", g1[!is.na(extreme_id), uniqueN(extreme_id)], "\n")

# Show event_80_95_id durations for this grid cell
g1_events <- g1[!is.na(event_80_95_id), .N, by = event_80_95_id]
cat("Event durations (event_80_95_id):\n")
print(g1_events[, .N, by = .(duration = N)][order(-N)])

# Show first 3 events in detail
first_events <- g1_events[order(event_80_95_id)][1:min(3, nrow(g1_events))]
for (eid in first_events$event_80_95_id) {
  cat("\n  Event", eid, "(duration =", first_events[event_80_95_id == eid, N], "):\n")
  ev <- g1[event_80_95_id == eid, .(date, value, std_value, event_80_95_id,
                                     extreme_id, prec_val)]
  print(ev)
}
cat("\n")

# ── 12. SPATIAL CHANGES FILE ──────────────────────────────────
cat("=== 12. SPATIAL_CHANGES.RDS ===\n")
sc_file <- paste0(PATH_OUTPUT, 'spatial_changes.rds')
if (file.exists(sc_file)) {
  sc <- readRDS(sc_file)
  cat("Dimensions:", nrow(sc), "x", ncol(sc), "\n")
  cat("Columns   :", paste(colnames(sc), collapse = ", "), "\n")
  cat("Variables :", paste(unique(sc$variable), collapse = ", "), "\n")
  cat("Periods   :", paste(unique(sc$period), collapse = ", "), "\n")
  cat("Ratio summary (all):\n")
  print(summary(sc$ratio))
  cat("NaN count :", sum(is.nan(sc$ratio)), "\n")
  cat("Inf count :", sum(is.infinite(sc$ratio)), "\n")
  rm(sc)
} else {
  cat("File not found\n")
}

cat("\n=== END OF DIAGNOSTICS ===\n")
sink()
rm(exeves, prec, evap_grid, exeves_prec); gc()
cat("Diagnostics written to: diagnostics.txt\n")
