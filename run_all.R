# Master script to run the complete ExEvEs analysis workflow
# Replicates imarkonis/ithaca/projects/exeves/stable/czechia/ up to 07c
# Adapted for Europe/Mediterranean using GLEAM evaporation + MSWEP precipitation

cat("=====================================================\n")
cat("  ExEvEs Analysis - Europe / Mediterranean Region    \n")
cat("  Evaporation (GLEAM) + Precipitation (MSWEP)       \n")
cat("=====================================================\n\n")

start_time <- Sys.time()

#---------------------------------------------------------------
# STEP 0: Initialize paths, constants, palettes
#---------------------------------------------------------------
cat("STEP 00: Initializing paths & constants...\n")
cat("-------------------------------------------\n")
source("00_initialize.R")

# Verify raw data exist
stopifnot("Evaporation NetCDF not found" = file.exists(EVAP_NC_FILE))
stopifnot("Precipitation NetCDF not found" = file.exists(PREC_NC_FILE))
cat("  Both raw NetCDF files found.\n\n")

#---------------------------------------------------------------
# STEP 1: Read & align GLEAM + MSWEP data
#---------------------------------------------------------------
cat("STEP 01: Reading & aligning NetCDF data...\n")
cat("-------------------------------------------\n")
cat("  (this may take 10-30 min for large grids)\n")
source("01_read_data.R")
cat("\n")

#---------------------------------------------------------------
# STEP 2: ExEvE identification (4 definitions)
#---------------------------------------------------------------
cat("STEP 02: Preprocessing – ExEvE identification...\n")
cat("-------------------------------------------\n")
cat("  (pentad standardisation + quantile regression)\n")
source("02_preprocessing.R")
cat("\n")

#---------------------------------------------------------------
# STEP 3: Temporal clustering (ACF analysis)
#---------------------------------------------------------------
cat("STEP 03: Temporal clustering...\n")
cat("-------------------------------------------\n")
source("03_temporal_clustering.R")
cat("\n")

#---------------------------------------------------------------
# STEP 4: Statistical properties of event definitions
#---------------------------------------------------------------
cat("STEP 04: Event statistical properties...\n")
cat("-------------------------------------------\n")
source("04_stat_properties.R")
cat("\n")

#---------------------------------------------------------------
# STEP 5a: Monthly changes (evap + prec)
#---------------------------------------------------------------
cat("STEP 05a: Monthly changes...\n")
cat("-------------------------------------------\n")
source("05a_changes_monthly.R")
cat("\n")

#---------------------------------------------------------------
# STEP 5b: Spatial changes between periods
#---------------------------------------------------------------
cat("STEP 05b: Spatial changes...\n")
cat("-------------------------------------------\n")
source("05b_changes_spatial.R")
cat("\n")

#---------------------------------------------------------------
# STEP 5c: Change visualisation (timeseries + polar + spatial)
#---------------------------------------------------------------
cat("STEP 05c: Change plots...\n")
cat("-------------------------------------------\n")
source("05c_changes_plots.R")
cat("\n")

#---------------------------------------------------------------
# STEP 7a: P-E implications (waffle plots)
#---------------------------------------------------------------
cat("STEP 07a: P-E implications...\n")
cat("-------------------------------------------\n")
source("07a_implications.R")
cat("\n")

#---------------------------------------------------------------
# STEP 7b: Monthly P-E implications
#---------------------------------------------------------------
cat("STEP 07b: Monthly P-E implications...\n")
cat("-------------------------------------------\n")
source("07b_implications_monthly.R")
cat("\n")

#---------------------------------------------------------------
# STEP 7c: Day of extremes within events
#---------------------------------------------------------------
cat("STEP 07c: Day of extremes...\n")
cat("-------------------------------------------\n")
source("07c_day_of_extremes.R")
cat("\n")

#---------------------------------------------------------------
# STEP 8a: Spatial maps (event frequency, duration, change)
#---------------------------------------------------------------
cat("STEP 8a: Spatial maps...\n")
cat("-------------------------------------------\n")
source("03_spatial_maps.R")
cat("\n")

#---------------------------------------------------------------
# STEP 8b: Temporal analysis (annual trends, seasonal patterns)
#---------------------------------------------------------------
cat("STEP 8b: Temporal analysis...\n")
cat("-------------------------------------------\n")
source("04_temporal_analysis.R")
cat("\n")

#---------------------------------------------------------------
# SUMMARY
#---------------------------------------------------------------
end_time <- Sys.time()
elapsed  <- difftime(end_time, start_time, units = "mins")

cat("\n=====================================================\n")
cat("               ANALYSIS COMPLETE!                    \n")
cat("=====================================================\n\n")
cat("Total time:", round(as.numeric(elapsed), 1), "minutes\n\n")

cat("Output locations:\n")
cat("  Data    :", PATH_OUTPUT_DATA, "\n")
cat("  Figures :", PATH_OUTPUT_FIGURES, "\n")
cat("  Tables  :", PATH_OUTPUT_TABLES, "\n\n")

cat("Key outputs:\n")
cat("  data/exeves_std_europe_med.rds       – main ExEvE dataset\n")
cat("  data/europe_med_prec_grid.rds        – aligned precipitation grid\n")
cat("  tables/definition_properties.csv     – event property table\n")
cat("  tables/definition_change.csv         – period change table\n")
cat("  figures/clustering.png               – ACF / temporal clustering\n")
cat("  figures/exeve_changes.png            – severity/intensity changes\n")
cat("  figures/implications.png             – P-E budget (all/ExEvE/non)\n")
cat("  figures/implications_monthly.png     – monthly P-E budget\n")
cat("  figures/wet_days.png                 – day-of-extreme analysis\n")
cat("  figures/map_event_*.png              – spatial maps (frequency, duration, change)\n")
cat("  figures/timeseries_*.png             – temporal trends (annual, duration)\n\n")

cat("To explore results interactively:\n")
cat("  exeves <- readRDS('data/exeves_std_europe_med.rds')\n")
cat("  prec   <- readRDS('data/europe_med_prec_grid.rds')\n")
