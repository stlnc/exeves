# Extreme Evaporation Events (ExEvEs) Analysis for Europe/Mediterranean

This workflow analyzes extreme evaporation events based on the methodology from the ithaca project (imarkonis/ithaca on GitHub). It identifies periods of anomalously high evaporation and examines their spatial and temporal patterns.

## Data Requirements

You need gridded evaporation data in NetCDF format:
- **Current file**: `gleam_e_mm_europe_med.nc`
- **Variable name**: Should be "E" (evaporation)
- **Dimensions**: lon, lat, time
- **Units**: mm/day or mm (will be automatically detected)

## Workflow Steps

### 0. Initialize (`00_initialize.R`)
Sets up directory structure and paths:
- `data/` - Processed data files
- `figures/` - All plots and visualizations
- `tables/` - Summary statistics and tables

**Run this first!**

```r
source("00_initialize.R")
```

### 1. Read and Preprocess Data (`01_read_preprocess_data.R`)
**This is the core script** - it:
- Reads the NetCDF evaporation file
- Creates pentad-based climatologies (5-day periods)
- Standardizes values (z-scores)
- Identifies extreme events using multiple definitions:
  - **Q80/Q95**: Events are periods > 80th percentile containing ≥1 day > 95th percentile
  - **Mean/Q95**: Events are periods > mean containing ≥1 day > 95th percentile  
  - **Mean/Q95***: Same as Mean/Q95 but with quantile regression
  - **Q80**: Simple threshold (all days > 80th percentile)

**Outputs:**
- `data/europe_med_evap_grid.rds` - Raw evaporation data
- `data/grid_europe_med.rds` - Grid cell coordinates
- `data/pentads_std_europe_med.rds` - Pentad statistics
- `data/exeves_std_europe_med.rds` - **Main output** with all event identifications

**Time**: ~10-30 minutes depending on data size

```r
source("01_read_preprocess_data.R")
```

### 2. Statistical Properties (`02_stat_properties.R`)
Calculates and compares different event definitions:
- Event frequency (events per grid cell)
- Event duration (mean days per event)
- Event severity (total evaporation per event)
- Event intensity (mean evaporation rate during events)
- Temporal changes (1981-2001 vs 2002-2022)
- Seasonal distribution

**Outputs:**
- `tables/europe_med_event_properties.csv` - Comparison table
- Console output with detailed statistics

```r
source("02_stat_properties.R")
```

### 3. Spatial Maps (`03_spatial_maps.R`)
Creates maps showing:
- Event frequency across the region
- Mean event duration
- Change in frequency between periods
- Ratio of change (2002-2022 / 1981-2001)

**Outputs:**
- `figures/map_event_frequency.png`
- `figures/map_event_duration.png`
- `figures/map_event_change.png`
- `figures/map_event_ratio.png`

```r
source("03_spatial_maps.R")
```

### 4. Temporal Analysis (`04_temporal_analysis.R`)
Analyzes time series and trends:
- Annual event frequency over time
- Monthly climatology
- Seasonal patterns
- Duration trends
- Statistical significance of trends

**Outputs:**
- `figures/timeseries_annual.png`
- `figures/climatology_monthly.png`
- `figures/seasonal_comparison.png`
- `figures/timeseries_duration.png`
- `tables/temporal_summary.csv`

```r
source("04_temporal_analysis.R")
```

## Running the Complete Workflow

You can run all scripts in sequence:

```r
# Run the complete workflow
source("00_initialize.R")
source("01_read_preprocess_data.R")
source("02_stat_properties.R")
source("03_spatial_maps.R")
source("04_temporal_analysis.R")
```

Or use the master script:

```r
source("run_all.R")
```

## Understanding the Results

### Event Definitions

The analysis uses the **Q80/Q95 definition** as the primary one (following the ithaca project):
- An event starts when evaporation > 80th percentile
- The event must contain at least one day with evaporation > 95th percentile
- The event ends when evaporation drops below the 80th percentile

This captures sustained periods of high evaporation (not just individual extreme days).

### Key Metrics

- **Frequency**: Number of events per grid cell over the study period
- **Duration**: How many days each event lasts
- **Severity**: Total evaporation during the event (mm)
- **Intensity**: Mean evaporation rate during the event (mm/day)

### Time Periods

The analysis splits data into two periods:
- **Period 1**: 1981-2001 (21 years)
- **Period 2**: 2002-2022 (21 years)

This allows detection of temporal changes in extreme event characteristics.

## Required R Packages

```r
install.packages(c(
  "data.table",   # Data manipulation
  "ncdf4",        # NetCDF file handling
  "lubridate",    # Date handling
  "quantreg",     # Quantile regression
  "ggplot2",      # Plotting
  "sf"            # Spatial data (optional, for boundaries)
))
```

## Troubleshooting

### "Cannot open file" error
- Check that `gleam_e_mm_europe_med.nc` is in the working directory
- Verify the filename matches exactly

### Memory issues
- The script processes data in chunks, but very large datasets may still cause issues
- Try closing other programs
- Consider processing smaller regions

### NetCDF variable name issues
- Open the file and check variable names:
  ```r
  library(ncdf4)
  nc <- nc_open("gleam_e_mm_europe_med.nc")
  names(nc$var)
  nc_close(nc)
  ```
- Modify script line `evap_array <- ncvar_get(nc, "E")` to match your variable name

### Time unit issues
- The script auto-detects "days since" or "hours since" formats
- If dates are wrong, check the NetCDF time units and modify the conversion logic

## Citation

This workflow is adapted from the ITHACA project:
- Repository: https://github.com/imarkonis/ithaca
- Methods based on extreme event identification in hydroclimatology

## Next Steps

After running this basic workflow, you can:

1. **Add additional variables**: precipitation, radiation, temperature to understand drivers
2. **Compute event drivers**: What meteorological conditions lead to extreme events?
3. **Regional analysis**: Split by climate zones or sub-regions
4. **Advanced statistics**: Trend attribution, circulation patterns
5. **Compare datasets**: Run with different evaporation products (GLEAM, ERA5, etc.)

## Contact & Support

For questions about:
- **This workflow**: Check the ithaca repository documentation
- **GLEAM data**: https://www.gleam.eu/
- **Methodology**: See publications from the ithaca project

---

**Last updated**: 2026-02-25
