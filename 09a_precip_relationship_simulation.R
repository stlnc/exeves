# Simulation-based uncertainty analysis of ExEvE-precipitation relationship
#
# Goal:
# 1) Test timing relationship using lag = (prec peak day - first extreme evap day)
# 2) Test mean precipitation anomaly during ExEvEs relative to matched non-ExEvE controls
# 3) Quantify uncertainty via Monte Carlo null simulations + bootstrap spread
#
# Outputs:
# - figures/precip_relationship_simulation.png
# - tables/precip_relationship_summary.csv
# - data/precip_relationship_null_bootstrap.rds

library(data.table)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(parallel)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")
load(paste0(PATH_OUTPUT_DATA, "grid_cell_n.Rdata"))

#===============================================================================
# SETTINGS
#===============================================================================
N_SIM <- 2000L
N_BOOT <- 1000L
MIN_DURATION <- 3L
DURATION_TOL <- 2L
MAX_EVENTS_INFERENCE <- suppressWarnings(as.integer(Sys.getenv("EXEVES_MAX_EVENTS_INFERENCE", unset = "200000")))

# Use scheduler allocation when available; override with EXEVES_N_CORES.
# Override with Sys.setenv(EXEVES_N_CORES = "<n>")
cores_env <- suppressWarnings(as.integer(Sys.getenv("EXEVES_N_CORES", unset = "")))
sched_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
if (is.na(sched_cores)) sched_cores <- suppressWarnings(as.integer(Sys.getenv("PBS_NP", unset = "")))
if (is.na(sched_cores)) sched_cores <- suppressWarnings(as.integer(Sys.getenv("NSLOTS", unset = "")))
if (is.na(sched_cores)) sched_cores <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "")))
detected <- parallel::detectCores(logical = TRUE)
if (is.na(detected) || detected < 1L) detected <- 1L
if (is.na(cores_env)) {
  if (is.na(sched_cores)) {
    N_CORES <- max(1L, detected - 1L)
  } else {
    N_CORES <- max(1L, min(sched_cores, detected))
  }
} else {
  N_CORES <- max(1L, min(cores_env, detected))
}

cat("Running 09a simulation analysis...\n")
cat("  N_SIM      :", N_SIM, "\n")
cat("  N_BOOT     :", N_BOOT, "\n")
cat("  Max events :", MAX_EVENTS_INFERENCE, "\n")
cat("  Cores used :", N_CORES, "/", detected, "\n")

set.seed(42)

#===============================================================================
# 1. LOAD + MERGE DATA
#===============================================================================
cat("Loading data...\n")
exeves <- readRDS(paste0(PATH_OUTPUT_DATA, "exeves_std_", region, ".rds"))
prec   <- readRDS(paste0(PATH_OUTPUT_DATA, region, "_prec_grid.rds"))

setkey(exeves, grid_id, date)
setkey(prec,   grid_id, date)
exeves[prec, prec := i.value, on = .(grid_id, date)]
setnames(exeves, "value", "evap")
rm(prec); gc()

#===============================================================================
# 2. EVENT METRICS
#===============================================================================
cat("Computing event-level metrics...\n")

evt <- exeves[!is.na(event_80_95_id)][order(grid_id, event_80_95_id, date)]
evt[, event_day := seq_len(.N), by = .(grid_id, event_80_95_id)]

# First day with extreme evaporation in each event
first_ext <- evt[!is.na(extreme_id), .(event_day_of_first_extreme = min(event_day)),
                 by = .(grid_id, event_80_95_id)]

# Day of precipitation maximum within event
peak_idx <- evt[, .I[which.max(prec)], by = .(grid_id, event_80_95_id)]$V1
prec_peak <- evt[peak_idx, .(grid_id, event_80_95_id, event_day_of_first_prec_max = event_day)]

# Event boundaries and attributes
events <- evt[, .(
  start_date = min(date),
  end_date   = max(date),
  duration   = .N,
  month      = month(min(date)),
  period     = period[1],
  mean_prec_event = mean(prec, na.rm = TRUE)
), by = .(grid_id, event_80_95_id)]

events <- merge(events, first_ext, by = c("grid_id", "event_80_95_id"), all.x = TRUE)
events <- merge(events, prec_peak, by = c("grid_id", "event_80_95_id"), all.x = TRUE)
events[, lag := event_day_of_first_prec_max - event_day_of_first_extreme]

# Keep analysable events
events <- events[duration >= MIN_DURATION & !is.na(event_day_of_first_extreme)]
cat("  Events retained:", nrow(events), "\n")

if (nrow(events) == 0L) {
  stop("No events remain after filtering. Reduce MIN_DURATION or check inputs.")
}

# For Monte Carlo/Bootstrap, cap event count to keep runtime tractable on HPC.
events_inf <- events
if (!is.na(MAX_EVENTS_INFERENCE) && MAX_EVENTS_INFERENCE > 0L && nrow(events_inf) > MAX_EVENTS_INFERENCE) {
  cat("  Downsampling events for inference to:", MAX_EVENTS_INFERENCE, "\n")
  frac <- MAX_EVENTS_INFERENCE / nrow(events_inf)
  events_inf <- events_inf[, .SD[sample.int(.N, max(1L, floor(.N * frac)))], by = .(period, month)]
  if (nrow(events_inf) > MAX_EVENTS_INFERENCE) {
    events_inf <- events_inf[sample.int(.N, MAX_EVENTS_INFERENCE)]
  }
  cat("  Events used for inference:", nrow(events_inf), "\n")
}

#===============================================================================
# 3. BUILD NON-EXEVE CONTROL POOL (by grid_id x month, contiguous runs)
#===============================================================================
cat("Building matched control pool...\n")
non_evt <- exeves[is.na(event_80_95_id), .(grid_id, date, prec)]
non_evt[, month := month(date)]
setorder(non_evt, grid_id, month, date)
non_evt[, day_num := as.integer(date)]
non_evt[, run_id := rleid(day_num - seq_len(.N)), by = .(grid_id, month)]

run_tbl <- non_evt[, .(
  len = .N,
  prec_vec = list(prec)
), by = .(grid_id, month, run_id)]

run_tbl[, key := paste(grid_id, month, sep = "_")]
run_pool <- split(run_tbl[, .(key, len, prec_vec)], by = "key", keep.by = FALSE)

rm(non_evt, run_tbl); gc()

sample_control_mean_prec <- function(grid_id, month, duration, pool, tol = 2L) {
  key <- paste(grid_id, month, sep = "_")
  runs <- pool[[key]]
  if (is.null(runs) || nrow(runs) == 0L) return(NA_real_)

  # Strict duration first, then fallback to duration-1, duration-2, ...
  for (d in seq(from = duration, to = max(1L, duration - tol), by = -1L)) {
    eligible <- which(runs$len >= d)
    if (length(eligible) == 0L) next

    i_run <- if (length(eligible) == 1L) eligible else sample(eligible, 1L)
    vec <- runs$prec_vec[[i_run]]
    max_start <- runs$len[i_run] - d + 1L
    if (max_start < 1L) next
    i_start <- if (max_start == 1L) 1L else sample.int(max_start, 1L)
    return(mean(vec[i_start:(i_start + d - 1L)], na.rm = TRUE))
  }

  NA_real_
}

#===============================================================================
# 4. OBSERVED EFFECTS
#===============================================================================
cat("Computing observed effects...\n")
obs_lag_median <- events_inf[!is.na(lag), median(lag, na.rm = TRUE)]

set.seed(123)
obs_ctrl <- mapply(
  FUN = sample_control_mean_prec,
  grid_id = events_inf$grid_id,
  month = events_inf$month,
  duration = events_inf$duration,
  MoreArgs = list(pool = run_pool, tol = DURATION_TOL)
)
obs_deltaP <- mean(events_inf$mean_prec_event - obs_ctrl, na.rm = TRUE)

cat("  Observed median lag:", round(obs_lag_median, 3), "\n")
cat("  Observed mean DeltaP:", round(obs_deltaP, 4), "\n")

#===============================================================================
# 5. PARALLEL SIMULATION + BOOTSTRAP
#===============================================================================
simulate_one <- function(seed, events_dt, pool, tol) {
  set.seed(seed)

  # Null lag: random precip-peak day within each event duration
  null_peak <- vapply(events_dt$duration, function(d) sample.int(d, 1L), integer(1))
  null_lag <- null_peak - events_dt$event_day_of_first_extreme
  null_lag_median <- median(null_lag, na.rm = TRUE)

  # Null DeltaP: event mean prec vs matched non-event random windows
  ctrl <- mapply(
    FUN = sample_control_mean_prec,
    grid_id = events_dt$grid_id,
    month = events_dt$month,
    duration = events_dt$duration,
    MoreArgs = list(pool = pool, tol = tol)
  )
  null_deltaP_mean <- mean(events_dt$mean_prec_event - ctrl, na.rm = TRUE)

  c(null_lag_median = null_lag_median, null_deltaP_mean = null_deltaP_mean)
}

bootstrap_one <- function(seed, events_dt, pool, tol) {
  set.seed(seed)
  idx <- sample.int(nrow(events_dt), nrow(events_dt), replace = TRUE)
  samp <- events_dt[idx]

  lag_med <- median(samp$lag, na.rm = TRUE)
  ctrl <- mapply(
    FUN = sample_control_mean_prec,
    grid_id = samp$grid_id,
    month = samp$month,
    duration = samp$duration,
    MoreArgs = list(pool = pool, tol = tol)
  )
  dP <- mean(samp$mean_prec_event - ctrl, na.rm = TRUE)

  c(boot_lag_median = lag_med, boot_deltaP_mean = dP)
}

run_in_parallel_with_progress <- function(seeds, worker_fun, label, chunk_size = NULL) {
  if (is.null(chunk_size)) {
    chunk_size <- max(20L, N_CORES * 5L)
  }

  chunks <- split(seeds, ceiling(seq_along(seeds) / chunk_size))
  n_chunks <- length(chunks)
  out_chunks <- vector("list", n_chunks)

  cat(label, "\n")
  pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  t0 <- Sys.time()

  for (i in seq_len(n_chunks)) {
    out_chunks[[i]] <- mclapply(chunks[[i]], worker_fun,
                                mc.cores = N_CORES, mc.set.seed = FALSE)
    setTxtProgressBar(pb, i)

    if (i %% 5L == 0L || i == n_chunks) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rate <- i / max(elapsed, 1e-6)
      eta <- (n_chunks - i) / max(rate, 1e-6)
      cat(sprintf("  %s: chunk %d/%d | elapsed %.1f min | ETA %.1f min\n",
                  label, i, n_chunks, elapsed / 60, eta / 60))
    }
  }

  close(pb)
  do.call(c, out_chunks)
}

sim_seeds <- 100000L + seq_len(N_SIM)
sim_out <- run_in_parallel_with_progress(
  sim_seeds,
  function(s) simulate_one(s, events_inf, run_pool, DURATION_TOL),
  label = "Monte Carlo"
)
null_dt <- as.data.table(do.call(rbind, sim_out))

boot_seeds <- 200000L + seq_len(N_BOOT)
boot_out <- run_in_parallel_with_progress(
  boot_seeds,
  function(s) bootstrap_one(s, events_inf, run_pool, DURATION_TOL),
  label = "Bootstrap"
)
boot_dt <- as.data.table(do.call(rbind, boot_out))

#===============================================================================
# 6. INFERENCE SUMMARY
#===============================================================================
# Two-sided p-values against null distributions
p_lag <- mean(abs(null_dt$null_lag_median) >= abs(obs_lag_median), na.rm = TRUE)
p_dP  <- mean(abs(null_dt$null_deltaP_mean) >= abs(obs_deltaP), na.rm = TRUE)

ci_lag <- quantile(boot_dt$boot_lag_median, probs = c(0.025, 0.975), na.rm = TRUE)
ci_dP  <- quantile(boot_dt$boot_deltaP_mean, probs = c(0.025, 0.975), na.rm = TRUE)

summary_dt <- data.table(
  metric = c("median_lag_days", "mean_deltaP_mm_day"),
  observed = c(obs_lag_median, obs_deltaP),
  null_mean = c(mean(null_dt$null_lag_median, na.rm = TRUE),
                mean(null_dt$null_deltaP_mean, na.rm = TRUE)),
  null_sd = c(sd(null_dt$null_lag_median, na.rm = TRUE),
              sd(null_dt$null_deltaP_mean, na.rm = TRUE)),
  p_value_two_sided = c(p_lag, p_dP),
  ci_2.5 = c(ci_lag[1], ci_dP[1]),
  ci_97.5 = c(ci_lag[2], ci_dP[2])
)

cat("\n--- Simulation summary ---\n")
print(summary_dt)

#===============================================================================
# 7. PLOTS (spread-focused)
#===============================================================================
cat("Creating plots...\n")

# Panel A: Null distribution of lag statistic
p1 <- ggplot(null_dt, aes(x = null_lag_median)) +
  geom_histogram(bins = 45, fill = PALETTES$subdued_prof[1], color = "white") +
  geom_vline(xintercept = obs_lag_median, color = PALETTES$subdued_prof[4], linewidth = 0.9) +
  xlab("Null median lag (prec peak day - extreme evap day)") +
  ylab("Simulation count") +
  labs(subtitle = paste0(
    "Observed=", round(obs_lag_median, 2),
    " | p=", signif(p_lag, 3)
  )) +
  theme_linedraw()

# Panel B: Null distribution of DeltaP statistic
p2 <- ggplot(null_dt, aes(x = null_deltaP_mean)) +
  geom_histogram(bins = 45, fill = PALETTES$subdued_prof[1], color = "white") +
  geom_vline(xintercept = obs_deltaP, color = PALETTES$subdued_prof[4], linewidth = 0.9) +
  xlab("Null mean DeltaP (event - matched control, mm/day)") +
  ylab("Simulation count") +
  labs(subtitle = paste0(
    "Observed=", round(obs_deltaP, 3),
    " | p=", signif(p_dP, 3)
  )) +
  theme_linedraw()

# Panel C: Bootstrap spread with 95% CI
boot_long <- rbind(
  data.table(metric = "Median lag (days)", value = boot_dt$boot_lag_median,
             observed = obs_lag_median, ci_lo = ci_lag[1], ci_hi = ci_lag[2]),
  data.table(metric = "Mean DeltaP (mm/day)", value = boot_dt$boot_deltaP_mean,
             observed = obs_deltaP, ci_lo = ci_dP[1], ci_hi = ci_dP[2])
)

p3 <- ggplot(boot_long, aes(x = value)) +
  geom_histogram(bins = 45, fill = PALETTES$subdued_prof[2], color = "white") +
  geom_vline(aes(xintercept = observed), color = PALETTES$subdued_prof[4], linewidth = 0.9) +
  geom_vline(aes(xintercept = ci_lo), linetype = 2, color = "gray35") +
  geom_vline(aes(xintercept = ci_hi), linetype = 2, color = "gray35") +
  facet_wrap(~ metric, scales = "free_x", ncol = 1) +
  xlab("Bootstrap statistic value") +
  ylab("Bootstrap count") +
  labs(subtitle = "Dashed lines = 95% bootstrap CI") +
  theme_linedraw()

cat("Saving figure...\n")
ggarrange(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"))
ggsave(paste0(PATH_OUTPUT_FIGURES, "precip_relationship_simulation.png"),
       width = 10, height = 14, dpi = 300)

#===============================================================================
# 8. SAVE OUTPUT TABLES
#===============================================================================
write.csv(summary_dt,
          paste0(PATH_OUTPUT_TABLES, "precip_relationship_summary.csv"),
          row.names = FALSE)

saveRDS(
  list(
    settings = list(N_SIM = N_SIM, N_BOOT = N_BOOT, N_CORES = N_CORES,
                    MIN_DURATION = MIN_DURATION, DURATION_TOL = DURATION_TOL,
                    MAX_EVENTS_INFERENCE = MAX_EVENTS_INFERENCE),
    observed = list(obs_lag_median = obs_lag_median, obs_deltaP = obs_deltaP),
    null = null_dt,
    bootstrap = boot_dt,
    summary = summary_dt
  ),
  paste0(PATH_OUTPUT_DATA, "precip_relationship_null_bootstrap.rds")
)

cat("Saved:", paste0(PATH_OUTPUT_FIGURES, "precip_relationship_simulation.png"), "\n")
cat("Saved:", paste0(PATH_OUTPUT_TABLES, "precip_relationship_summary.csv"), "\n")
cat("Saved:", paste0(PATH_OUTPUT_DATA, "precip_relationship_null_bootstrap.rds"), "\n")

rm(exeves, evt, events, run_pool, null_dt, boot_dt, boot_long, summary_dt); gc()
cat("09a simulation analysis complete.\n")
