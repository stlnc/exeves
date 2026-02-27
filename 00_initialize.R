# Creates project paths for Europe/Mediterranean exeves analysis

# Note: Adjust these paths based on your local setup
# This assumes you're working in the exeves directory

## Paths
PATH_SAVE <- getwd()  # Current working directory
PATH_OUTPUT <- paste0(PATH_SAVE, "/")
PATH_OUTPUT_DATA <- paste0(PATH_OUTPUT, "data/")
PATH_OUTPUT_RAW <- paste0(PATH_OUTPUT_DATA, "raw/")
PATH_OUTPUT_FIGURES <- paste0(PATH_OUTPUT, "figures/")
PATH_OUTPUT_TABLES <- paste0(PATH_OUTPUT, "tables/")

dir.create(PATH_OUTPUT_DATA, showWarnings = FALSE)
dir.create(PATH_OUTPUT_RAW, showWarnings = FALSE)
dir.create(PATH_OUTPUT_FIGURES, showWarnings = FALSE)
dir.create(PATH_OUTPUT_TABLES, showWarnings = FALSE)

save(PATH_OUTPUT,
     PATH_OUTPUT_RAW,
     PATH_OUTPUT_DATA,
     PATH_OUTPUT_FIGURES,
     PATH_OUTPUT_TABLES,
     file = paste0(PATH_OUTPUT, "paths.Rdata"))

cat("Paths initialized successfully!\n")
cat("Data directory:", PATH_OUTPUT_DATA, "\n")
cat("Figures directory:", PATH_OUTPUT_FIGURES, "\n")
