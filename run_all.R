# Master script to run complete ExEvEs analysis workflow
# Runs all analysis steps in sequence

cat("===========================================\n")
cat("  ExEvEs Analysis - Europe/Mediterranean  \n")
cat("===========================================\n\n")

start_time <- Sys.time()

# Check if data file exists
if (!file.exists("gleam_e_mm_europe_med.nc")) {
  stop("ERROR: Cannot find 'gleam_e_mm_europe_med.nc' in current directory!\n",
       "Please make sure the NetCDF file is in: ", getwd())
}

cat("Data file found: gleam_e_mm_europe_med.nc\n\n")

# Step 0: Initialize
cat("STEP 0: Initializing paths...\n")
cat("-------------------------------------------\n")
source("00_initialize.R")
cat("\n")

# Step 1: Read and preprocess
cat("\nSTEP 1: Reading and preprocessing data...\n")
cat("-------------------------------------------\n")
cat("This may take 10-30 minutes...\n")
source("01_read_preprocess_data.R")
cat("\n")

# Step 2: Statistical properties
cat("\nSTEP 2: Computing statistical properties...\n")
cat("-------------------------------------------\n")
source("02_stat_properties.R")
cat("\n")

# Step 3: Spatial maps
cat("\nSTEP 3: Creating spatial maps...\n")
cat("-------------------------------------------\n")
source("03_spatial_maps.R")
cat("\n")

# Step 4: Temporal analysis
cat("\nSTEP 4: Analyzing temporal patterns...\n")
cat("-------------------------------------------\n")
source("04_temporal_analysis.R")
cat("\n")

# Summary
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n===========================================\n")
cat("           ANALYSIS COMPLETE!              \n")
cat("===========================================\n\n")
cat("Total time:", round(elapsed, 1), "minutes\n\n")

cat("Output locations:\n")
cat("  - Data: data/\n")
cat("  - Figures: figures/\n")
cat("  - Tables: tables/\n\n")

cat("Key outputs:\n")
cat("  - Main data: data/exeves_std_europe_med.rds\n")
cat("  - Event stats: tables/europe_med_event_properties.csv\n")
cat("  - All figures: figures/*.png\n\n")

cat("To explore results:\n")
cat("  - View PNG files in figures/ folder\n")
cat("  - Open CSV files in tables/ folder\n")
cat("  - Load RDS files in R for further analysis:\n")
cat("    exeves <- readRDS('data/exeves_std_europe_med.rds')\n\n")
