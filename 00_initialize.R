# Creates project paths and defines constants for Europe/Mediterranean ExEvEs analysis
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/

#===============================================================================
# PATHS
#===============================================================================
PATH_SAVE <- getwd()
PATH_OUTPUT <- paste0(PATH_SAVE, "/")
PATH_OUTPUT_DATA <- paste0(PATH_OUTPUT, "data/")
PATH_OUTPUT_RAW <- paste0(PATH_OUTPUT_DATA, "raw/")
PATH_OUTPUT_FIGURES <- paste0(PATH_OUTPUT, "figures/")
PATH_OUTPUT_TABLES <- paste0(PATH_OUTPUT, "tables/")

dir.create(PATH_OUTPUT_DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_OUTPUT_RAW, showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_OUTPUT_FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(PATH_OUTPUT_TABLES, showWarnings = FALSE, recursive = TRUE)

save(PATH_OUTPUT, PATH_OUTPUT_RAW, PATH_OUTPUT_DATA, PATH_OUTPUT_FIGURES, PATH_OUTPUT_TABLES,
     file = paste0(PATH_OUTPUT, "paths.Rdata"))

#===============================================================================
# REGION
#===============================================================================
region <- 'europe_med'

#===============================================================================
# TIME PERIODS
# GLEAM evap: 1980-01 to 2024-12
# MSWEP prec: 1979-01 to 2020-12
# Common analysis period: 1981-01-01 to 2020-12-31
#===============================================================================
START_PERIOD_1 <- as.Date("1981-01-01")
END_PERIOD_1   <- as.Date("2001-12-31")
END_PERIOD_2   <- as.Date("2020-12-31")
PERIOD_LENGTH   <- round(as.numeric((END_PERIOD_2 - START_PERIOD_1) / 365.25), 0)

#===============================================================================
# THRESHOLDS
#===============================================================================
EXTREMES_THRES  <- 0.95
LOW_THRES       <- 0.80
SUB_PERIOD_YEARS <- 0.5 * PERIOD_LENGTH
DAYS_IN_YEAR    <- 365.25
SEC_IN_DAY      <- 86400

#===============================================================================
# COLOUR PALETTES
#===============================================================================
colset_subdued_prof   <- c("#90AFC5", "#336B87", "#2A3132", "#763626")
SUBDUED_PROF_PALETTE  <- colset_subdued_prof
agu_palette           <- c('#00324A', '#005294', '#058ECD', '#FFFFFF')
WATER_CYCLE_CHANGE_PALETTE <- c('steelblue3', 'darkgreen', 'darkred', 'darkorange')
# Order: Wetter-Accelerated, Wetter-Decelerated, Drier-Accelerated, Drier-Decelerated
PALETTES <- list(
  subdued_prof       = colset_subdued_prof,
  water_cycle_change = WATER_CYCLE_CHANGE_PALETTE,
  agu                = agu_palette
)

#===============================================================================
# RAW DATA FILES
#===============================================================================
EVAP_NC_FILE <- paste0(PATH_OUTPUT_RAW, "gleam_e_mm_med_198001_202412_025_daily.nc")
PREC_NC_FILE <- paste0(PATH_OUTPUT_RAW, "mswep-v2-8_tp_mm_med_197901_202012_025_daily.nc")

cat("Paths initialized successfully!\n")
cat("  Data     :", PATH_OUTPUT_DATA, "\n")
cat("  Figures  :", PATH_OUTPUT_FIGURES, "\n")
cat("  Tables   :", PATH_OUTPUT_TABLES, "\n")
cat("  Region   :", region, "\n")
cat("  Period   :", as.character(START_PERIOD_1), "to", as.character(END_PERIOD_2), "\n")
cat("  Split at :", as.character(END_PERIOD_1), "\n")
