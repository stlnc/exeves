# Onset/termination analysis: how drivers differ at start and end of events
# Adapted from imarkonis/ithaca/projects/exeves/stable/czechia/06c_drivers_onset_termination.R
# Uses precipitation as the available co-variable (radiation/temp not available)

library(data.table)
library(ggplot2)
library(lubridate)

# Load paths and constants
load("paths.Rdata")
source("00_initialize.R")

cat("Loading drivers data...\n")
exeves_drivers <- readRDS(paste0(PATH_OUTPUT_DATA, region, '_exeves_drivers.rds'))

#===============================================================================
# ONSET / TERMINATION DATES
#===============================================================================
onset_date <- exeves_drivers[event_day == 1 & conditions == 'ExEvE',
                              .(date = date), grid_id]
termination_date <- exeves_drivers[event_day == event_duration,
                                    .(date = date + 1), grid_id]
termination_date <- termination_date[date <= END_PERIOD_2]

# Drivers for analysis
exeves_drivers <- exeves_drivers[, .(grid_id, date, conditions, evap, prec)]

# Non-ExEvE baseline means by month
non_exeves_values <- melt(exeves_drivers[conditions == 'non-ExEvE'],
                          id.vars = c('grid_id', 'date', 'conditions'))
non_exeves_means <- non_exeves_values[, mean(value),
                                       .(month = month(date), conditions, variable)]

# Onset values
onset_values <- merge(onset_date, exeves_drivers, by = c('grid_id', 'date'), all.x = TRUE)
onset_values[, month := month(date)][, date := NULL]
onset_values <- melt(onset_values, id.vars = c('grid_id', 'month', 'conditions'))
onset_means <- onset_values[, mean(value), .(month, conditions, variable)]
levels(onset_means$conditions)[1] <- "Onset"

# Termination values
termination_values <- merge(termination_date, exeves_drivers,
                            by = c('grid_id', 'date'), all.x = TRUE)
termination_values[, month := month(date)][, date := NULL]
termination_values <- melt(termination_values, id.vars = c('grid_id', 'month', 'conditions'))
termination_means <- termination_values[, mean(value), .(month, conditions, variable)]
levels(termination_means$conditions)[2] <- "Termination"

# Combine
onset_termination <- rbind(onset_means, termination_means, non_exeves_means)
names(onset_termination)[4] <- "value"
onset_termination$conditions <- factor(onset_termination$conditions,
                                        levels = c("Onset", "Termination", "non-ExEvE"))
onset_termination$variable <- factor(onset_termination$variable,
                                      levels = c("evap", "prec"),
                                      labels = c("'Evaporation (mm/day)'",
                                                 "'Precipitation (mm/day)'"))

#===============================================================================
# PLOT
#===============================================================================
cat("Creating onset/termination plot...\n")

ggplot(onset_termination) +
  geom_line(aes(y = value, x = factor(month), col = conditions, group = conditions)) +
  geom_point(aes(y = value, x = factor(month), col = conditions, group = conditions)) +
  geom_line(data = onset_termination[conditions != "non-ExEvE"],
            aes(y = value, x = factor(month), group = month),
            col = SUBDUED_PROF_PALETTE[2], lty = 3) +
  facet_wrap(~variable, scales = "free", strip.position = "left",
             labeller = label_parsed) +
  labs(colour = "Conditions") +
  xlab("Month") + ylab("") +
  scale_color_manual(values = c(SUBDUED_PROF_PALETTE[c(2, 1)], "grey50")) +
  theme_linedraw() +
  theme(panel.grid.minor = element_line(colour = "grey60"),
        panel.grid.major = element_line(colour = "grey60"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(colour = 'black'))
ggsave(paste0(PATH_OUTPUT_FIGURES, "onset_termination.png"), width = 8, height = 4)

cat("Onset/termination analysis complete.\n")
rm(exeves_drivers); gc()
