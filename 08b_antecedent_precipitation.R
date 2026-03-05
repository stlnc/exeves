# Antecedent and post-event precipitation analysis
# Evaluates cumulative precipitation in the 7 days before each ExEvE starts
# and the 7 days after each ExEvE ends, to assess:
#   1. Whether pre-event moisture primes extreme evaporation
#   2. Whether evaporated water contributes to subsequent precipitation
#
# Panel A: Composite daily precipitation & evaporation ±7 days around ExEvE start
# Panel B: 7-day pre-event cumulative precipitation by month and period
# Panel C: 7-day post-event cumulative precipitation by month and period
#
# Output: figures/antecedent_precipitation.png

library(data.table)
library(ggplot2)
library(ggpubr)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, 'grid_cell_n.Rdata'))

WINDOW <- 7L  # days before start / after end

cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, 'exeves_std_', region, '.rds'))
prec   <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_prec_grid.rds'))

# Build a full daily lookup with both evap and prec
setnames(exeves, "value", "evap")
setnames(prec,   "value", "prec")
setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.prec, on = .(grid_id, date)]
rm(prec); gc()

#===============================================================================
# Identify event boundaries
#===============================================================================
cat("Identifying event boundaries...\n")

evt_rows <- exeves[!is.na(event_80_95_id)]
events <- evt_rows[, .(start_date = min(date),
                        end_date   = max(date),
                        period     = period[1],
                        month      = month(min(date)),
                        duration   = .N),
                    by = .(grid_id, event_80_95_id)]
rm(evt_rows); gc()
cat("  Total events:", nrow(events), "\n")

#===============================================================================
# Panel A: Composite daily profile around event start
#
# Memory-efficient approach: iterate through offsets, join + aggregate
# per-offset, only keeping the summary table.
#===============================================================================
cat("Panel A: Building composite daily profile...\n")

offsets_comp <- seq(-WINDOW, WINDOW)
comp_list <- vector("list", length(offsets_comp))

for (i in seq_along(offsets_comp)) {
  off <- offsets_comp[i]
  events[, lookup_date := start_date + off]
  events[exeves, `:=`(d_prec = i.prec, d_evap = i.evap),
         on = .(grid_id, lookup_date == date)]
  comp_list[[i]] <- events[, .(mean_prec = mean(d_prec, na.rm = TRUE),
                                mean_evap = mean(d_evap, na.rm = TRUE),
                                offset    = off),
                            by = .(month, period)]
  events[, c("d_prec", "d_evap", "lookup_date") := NULL]
}

composite_mean <- rbindlist(comp_list)
rm(comp_list); gc()

composite_mean[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
composite_mean[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]

comp_long <- melt(composite_mean,
                  id.vars = c("offset", "month", "period", "Period"),
                  measure.vars = c("mean_prec", "mean_evap"),
                  variable.name = "Variable", value.name = "value")
comp_long[, Variable := fifelse(Variable == "mean_prec", "Precipitation", "Evaporation")]

gg_composite <- ggplot(comp_long) +
  geom_vline(xintercept = 0, linetype = 2, col = 'red', alpha = 0.6) +
  geom_hline(yintercept = 0, col = 'black') +
  geom_line(aes(x = offset, y = value, col = Variable, linetype = Period),
            linewidth = 0.6) +
  facet_wrap(~ month, scales = 'free_y', ncol = 4) +
  scale_color_manual(values = c("Precipitation" = PALETTES$subdued_prof[2],
                                 "Evaporation"   = PALETTES$subdued_prof[4])) +
  scale_linetype_manual(values = c("Up to 2001" = 1, "After 2001" = 2)) +
  xlab("Days relative to ExEvE start") +
  ylab("Mean flux (mm/day)") +
  labs(subtitle = "Red dashed = ExEvE start day; composite across all event durations") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"),
        strip.background = element_rect(fill = 'grey30'))

rm(composite_mean, comp_long); gc()

#===============================================================================
# Panel B: 7-day cumulative precipitation BEFORE event start
#
# Same offset-loop approach: accumulate into events$prec_pre
#===============================================================================
cat("Panel B: Pre-event cumulative precipitation...\n")

events[, prec_pre := 0]
for (off in seq(-WINDOW, -1L)) {
  events[, lookup_date := start_date + off]
  events[exeves, daily_prec := i.prec, on = .(grid_id, lookup_date == date)]
  events[is.na(daily_prec), daily_prec := 0]
  events[, prec_pre := prec_pre + daily_prec]
  events[, c("daily_prec", "lookup_date") := NULL]
}

pre_cum <- events[, .(grid_id, event_80_95_id, period, month, prec_pre)]
pre_cum[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
pre_cum[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]

gg_pre <- ggplot(pre_cum) +
  geom_boxplot(aes(x = factor(month), y = prec_pre, fill = Period),
               outlier.size = 0.3, outlier.alpha = 0.15) +
  scale_fill_manual(values = c("Up to 2001" = PALETTES$subdued_prof[2],
                                "After 2001" = PALETTES$subdued_prof[4])) +
  xlab("Month") +
  ylab("Cumulative precipitation (mm)") +
  labs(subtitle = "Cumulative precipitation in the 7 days before ExEvE start") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"))

#===============================================================================
# Panel C: 7-day cumulative precipitation AFTER event end
#===============================================================================
cat("Panel C: Post-event cumulative precipitation...\n")

events[, prec_post := 0]
for (off in seq(1L, WINDOW)) {
  events[, lookup_date := end_date + off]
  events[exeves, daily_prec := i.prec, on = .(grid_id, lookup_date == date)]
  events[is.na(daily_prec), daily_prec := 0]
  events[, prec_post := prec_post + daily_prec]
  events[, c("daily_prec", "lookup_date") := NULL]
}

post_cum <- events[, .(grid_id, event_80_95_id, period, month, prec_post)]
post_cum[, Period := fifelse(period == "up_to_2001", "Up to 2001", "After 2001")]
post_cum[, Period := factor(Period, levels = c("Up to 2001", "After 2001"))]

gg_post <- ggplot(post_cum) +
  geom_boxplot(aes(x = factor(month), y = prec_post, fill = Period),
               outlier.size = 0.3, outlier.alpha = 0.15) +
  scale_fill_manual(values = c("Up to 2001" = PALETTES$subdued_prof[2],
                                "After 2001" = PALETTES$subdued_prof[4])) +
  xlab("Month") +
  ylab("Cumulative precipitation (mm)") +
  labs(subtitle = "Cumulative precipitation in the 7 days after ExEvE end") +
  theme_linedraw() +
  theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"))

#===============================================================================
# Print summary statistics
#===============================================================================
cat("\n--- Summary statistics ---\n")
cat("Pre-event cumulative precip (mm):\n")
print(pre_cum[, .(median = median(prec_pre), mean = mean(prec_pre),
                   q25 = quantile(prec_pre, 0.25), q75 = quantile(prec_pre, 0.75)),
               by = Period])
cat("\nPost-event cumulative precip (mm):\n")
print(post_cum[, .(median = median(prec_post), mean = mean(prec_post),
                    q25 = quantile(prec_post, 0.25), q75 = quantile(prec_post, 0.75)),
                by = Period])

# Difference (post - pre) per event
cum_both <- merge(pre_cum[, .(grid_id, event_80_95_id, prec_pre, Period)],
                  post_cum[, .(grid_id, event_80_95_id, prec_post)],
                  by = c("grid_id", "event_80_95_id"))
cum_both[, delta := prec_post - prec_pre]
cat("\nPost minus Pre (mm) — positive = more precip AFTER than BEFORE:\n")
print(cum_both[, .(median_delta = median(delta), mean_delta = mean(delta),
                    frac_positive = mean(delta > 0)),
                by = Period])

#===============================================================================
# Combine & save
#===============================================================================
cat("\nSaving figure...\n")
ggarrange(gg_composite, gg_pre, gg_post,
          ncol = 1, labels = c("A", "B", "C"),
          legend = 'bottom')
ggsave(paste0(PATH_OUTPUT_FIGURES, "antecedent_precipitation.png"),
       width = 10, height = 14, dpi = 300)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "antecedent_precipitation.png"), "\n")
rm(exeves, events, pre_cum, post_cum, cum_both); gc()
