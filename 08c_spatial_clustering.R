# Spatial clustering of ExEvEs
# For each date, identifies spatially connected clusters of grid cells
# that are simultaneously experiencing an ExEvE (event_80_95_id ≠ NA).
# Queen contiguity on the 0.25° grid defines adjacency.
#
# Panel A: Distribution of cluster sizes (# grid cells), by period
# Panel B: Mean daily cluster size by month, by period (boxplot)
# Panel C: Map of large-cluster participation frequency
#
# Output: figures/spatial_clustering.png

library(data.table)
library(ggplot2)
library(ggpubr)
library(igraph)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

#===============================================================================
# 1. LOAD DATA
#===============================================================================
cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
grid   <- readRDS(paste0(PATH_OUTPUT_DATA, 'grid_', region, '.rds'))

cat("  Grid cells:", nrow(grid), "\n")
cat("  Grid resolution: 0.25°\n")

#===============================================================================
# 2. PRECOMPUTE QUEEN-ADJACENCY EDGE LIST
#
# Two cells are queen-adjacent if |Δlon| ≤ 0.25 and |Δlat| ≤ 0.25
# (and not the same cell).  On a regular 0.25° grid this gives up
# to 8 neighbours per cell.
#===============================================================================
cat("Building queen-adjacency edge list...\n")

# Round to avoid floating-point near-misses
grid[, `:=`(lon_r = round(lon, 3), lat_r = round(lat, 3))]

# Self-join with inequality: CJ approach would be huge.
# Instead, for each cell enumerate its 8 potential neighbours and look them up.
offsets <- CJ(dlon = c(-0.25, 0, 0.25), dlat = c(-0.25, 0, 0.25))
offsets <- offsets[!(dlon == 0 & dlat == 0)]  # drop self

setkey(grid, lon_r, lat_r)
edges_list <- vector("list", nrow(offsets))

for (k in seq_len(nrow(offsets))) {
  tmp <- copy(grid[, .(grid_id, lon_r, lat_r)])
  tmp[, `:=`(lon_r = round(lon_r + offsets$dlon[k], 3),
             lat_r = round(lat_r + offsets$dlat[k], 3))]
  setkey(tmp, lon_r, lat_r)
  # Join: find grid cells whose shifted position matches another cell
  matched <- grid[tmp, nomatch = 0, on = .(lon_r, lat_r)]
  edges_list[[k]] <- matched[, .(from = i.grid_id, to = grid_id)]
}

edge_dt <- rbindlist(edges_list)
# Remove duplicates (a-b and b-a)
edge_dt[, `:=`(lo = pmin(from, to), hi = pmax(from, to))]
edge_dt <- unique(edge_dt[, .(from = lo, to = hi)])
cat("  Queen-adjacency edges:", nrow(edge_dt), "\n")

# Pre-build the full graph once (we'll take induced subgraphs per day)
g_full <- graph_from_data_frame(edge_dt, directed = FALSE,
                                 vertices = data.frame(name = grid$grid_id))

#===============================================================================
# 3. DAILY SPATIAL CLUSTERING
#
# For each date with at least one active ExEvE cell, find connected
# components in the induced subgraph of active cells.
#===============================================================================
cat("Running daily spatial clustering...\n")

# Active cell-dates: any cell currently inside an ExEvE
active <- exeves[!is.na(event_80_95_id), .(grid_id, date, period)]
active <- unique(active[, .(grid_id, date, period)])

# Group by date
active_by_date <- active[, .(grid_ids = list(grid_id),
                              period   = period[1]),
                          by = date]
n_dates <- nrow(active_by_date)
cat("  Dates with active ExEvEs:", n_dates, "\n")

# Pre-allocate result list
cluster_results <- vector("list", n_dates)

# Map grid_id → vertex index in g_full (igraph uses 1-based index)
vnames <- V(g_full)$name
id_to_idx <- setNames(seq_along(vnames), vnames)

t0 <- Sys.time()
for (i in seq_len(n_dates)) {
  gids <- active_by_date$grid_ids[[i]]
  n_active <- length(gids)

  if (n_active == 1L) {
    # Single cell = single cluster of size 1
    cluster_results[[i]] <- data.table(
      date       = active_by_date$date[i],
      period     = active_by_date$period[i],
      cluster_id = 1L,
      size       = 1L,
      n_active   = 1L,
      n_clusters = 1L
    )
    next
  }

  # Induced subgraph on active cells
  vidx <- id_to_idx[as.character(gids)]
  sg <- induced_subgraph(g_full, vidx)
  comp <- components(sg)

  sizes <- comp$csize
  n_cl  <- comp$no

  cluster_results[[i]] <- data.table(
    date       = active_by_date$date[i],
    period     = active_by_date$period[i],
    cluster_id = seq_len(n_cl),
    size       = sizes,
    n_active   = n_active,
    n_clusters = n_cl
  )

  if (i %% 2000 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  %d / %d dates  (%.0f s)\n", i, n_dates, elapsed))
  }
}

clusters <- rbindlist(cluster_results)
rm(cluster_results, active_by_date); gc()

cat("  Total cluster-date records:", nrow(clusters), "\n")
cat("  Median cluster size:", median(clusters$size), "\n")
cat("  Max cluster size:", max(clusters$size), "\n")

# Label periods for display
clusters[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
clusters[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]
clusters[, month := month(date)]

#===============================================================================
# 4. PER-CELL LARGE-CLUSTER PARTICIPATION
#
# For each cell, count how many days it belongs to a cluster with ≥ LARGE_THRESHOLD cells.
# This requires mapping cluster membership back to individual cells.
#===============================================================================
cat("Computing per-cell cluster participation...\n")
LARGE_THRESHOLD <- 10L  # cluster size to qualify as "large"

# Re-run components to get membership (cell → cluster_id mapping)
# More memory-efficient: iterate in batches
active2 <- exeves[!is.na(event_80_95_id), .(grid_id, date, period)]
active2 <- unique(active2[, .(grid_id, date)])
active2_by_date <- active2[, .(grid_ids = list(grid_id)), by = date]

cell_large_count <- data.table(grid_id = grid$grid_id, large_days = 0L)
setkey(cell_large_count, grid_id)

cat("  Mapping cell → cluster membership for large-cluster counting...\n")
t0 <- Sys.time()
for (i in seq_len(nrow(active2_by_date))) {
  gids <- active2_by_date$grid_ids[[i]]
  n_active <- length(gids)

  if (n_active < LARGE_THRESHOLD) next  # can't form a large cluster

  vidx <- id_to_idx[as.character(gids)]
  sg <- induced_subgraph(g_full, vidx)
  comp <- components(sg)

  # Which cells belong to large clusters?
  large_members <- which(comp$csize[comp$membership] >= LARGE_THRESHOLD)
  if (length(large_members) == 0L) next

  large_gids <- gids[large_members]
  cell_large_count[J(large_gids), large_days := large_days + 1L]

  if (i %% 2000 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  %d / %d dates  (%.0f s)\n", i, nrow(active2_by_date), elapsed))
  }
}

cell_large <- merge(grid, cell_large_count, by = "grid_id")
rm(active2, active2_by_date, cell_large_count); gc()

#===============================================================================
# 5. SUMMARY STATISTICS
#===============================================================================
cat("\n--- Summary statistics ---\n")
cat("Cluster size by period:\n")
print(clusters[, .(median_size = median(size),
                    mean_size   = mean(size),
                    q75         = quantile(size, 0.75),
                    q95         = quantile(size, 0.95),
                    max_size    = max(size),
                    n_clusters  = .N),
               by = Period])

cat("\nFraction of clusters with size >= 10:\n")
print(clusters[, .(frac_large = mean(size >= LARGE_THRESHOLD)), by = Period])

cat("\nTop 10 largest cluster-days:\n")
print(clusters[order(-size)][1:min(10, .N), .(date, Period, size, n_active, n_clusters)])

#===============================================================================
# 6. PLOTS
#===============================================================================
cat("Creating plots...\n")

# --- Panel A: ECDF of cluster sizes, by period ---
gg_size_dist <- ggplot(clusters) +
  stat_ecdf(aes(x = size, col = Period), linewidth = 0.7) +
  scale_color_manual(values = c("Up to 2001" = PALETTES$subdued_prof[2],
                                 "After 2001" = PALETTES$subdued_prof[4])) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  xlab("Cluster size (grid cells)") +
  ylab("Cumulative fraction") +
  labs(subtitle = "ECDF of daily spatial cluster sizes (queen contiguity)") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"))

# --- Panel B: Monthly mean cluster size, by period (boxplot) ---
# Summarise per date first: mean cluster size for that day
daily_summary <- clusters[, .(mean_size  = mean(size),
                               max_size   = max(size),
                               n_clusters = n_clusters[1]),
                           by = .(date, month, Period)]

gg_monthly <- ggplot(daily_summary) +
  geom_boxplot(aes(x = factor(month), y = mean_size, fill = Period),
               outlier.size = 0.3, outlier.alpha = 0.15) +
  scale_fill_manual(values = c("Up to 2001" = PALETTES$subdued_prof[2],
                                "After 2001" = PALETTES$subdued_prof[4])) +
  xlab("Month") +
  ylab("Mean cluster size (grid cells)") +
  labs(subtitle = "Daily mean spatial cluster size by month") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"))

# --- Panel C: Map of large-cluster participation frequency ---
gg_map <- ggplot(cell_large[large_days > 0]) +
  geom_tile(aes(x = lon, y = lat, fill = large_days)) +
  scale_fill_viridis_c(option = "inferno", name = "Days",
                        trans = "log1p",
                        breaks = c(1, 5, 20, 50, 100, 200)) +
  coord_quickmap() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(subtitle = paste0("Days in a large cluster (≥", LARGE_THRESHOLD, " connected cells)")) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        legend.position = "bottom",
        panel.grid = element_line(color = "gray90"))

#===============================================================================
# 7. COMBINE & SAVE
#===============================================================================
cat("Saving figure...\n")
ggarrange(gg_size_dist, gg_monthly, gg_map,
          ncol = 1, labels = c("A", "B", "C"),
          heights = c(1, 1, 1.3),
          legend = 'bottom')
ggsave(paste0(PATH_OUTPUT_FIGURES, "spatial_clustering.png"),
       width = 10, height = 14, dpi = 300)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "spatial_clustering.png"), "\n")

# Also save cluster data for potential further analysis
saveRDS(clusters, paste0(PATH_OUTPUT_DATA, "exeve_spatial_clusters.rds"))
cat("Saved cluster data:", paste0(PATH_OUTPUT_DATA, "exeve_spatial_clusters.rds"), "\n")

rm(exeves, grid, clusters, daily_summary, cell_large, g_full, edge_dt); gc()
