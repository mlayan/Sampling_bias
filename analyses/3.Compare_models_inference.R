#############################################################
##                  ANALYSIS OF DTA RUNS
##
## Creation Date: 13-05-2020
## Last Update: -05-2020
## Maylis Layan
#############################################################

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)

setwd("/mnt/gaia/MMMI_Rage/")
source('R_Functions/plot_results.R')

runType = 'mascot_v8'

#############################################################
## Load data
#############################################################
## Simulation data
m1 = read.table('HKY_M1/analyses/migration_events_sim1to10.txt', stringsAsFactors = FALSE, header = TRUE, sep = "\t")
m2 = read.table('HKY_M2/analyses/migration_events_sim1to10.txt', stringsAsFactors = FALSE, header = TRUE, sep = "\t")
m = bind_rows(m1, m2) %>%
  mutate(protocol = gsub("_[510]+$", "", protocol)) %>%
  filter(protocol != 'stratified') %>%
  mutate(protocol = factor(protocol, c(bias, surveillance))) %>%
  mutate(parameter = recode(parameter,
                            epidemicDuration = "treeHeight"))
## BEAST data
data = read.table(paste0("2.Figures/", runType, "/2.selected_data_", runType, ".txt"),
                  stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)),
         model = gsub("_.*$", "", model))

## Selected runs
selected = read.table(paste0("2.Figures/", runType, "/2.selected_runs_", runType, ".txt"),
                      stringsAsFactors = FALSE,
                      header = TRUE, sep = "\t")

#############################################################
##  Compare forwards/backwards migration rates with or
##  withouth BSSVS
#############################################################
## Migration rates with linear regression
spatialRates <- m %>%
  filter(parameter %in% c("locationFrequency", "migrationRate")) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(source = gsub("_[a-zA-Z]*$", "", parameterChar),
         dest = gsub("^[a-zA-Z]*_", "", parameterChar))

migrationRates <- data %>%
  # filter((grepl("^bssvs_", parameter) & model == "mascot"), BF >= 3) %>%
  # mutate(direction = gsub("^bssvs_|_[A-Z][A-Za-z_]*$", "", parameter)) %>%
  filter((grepl("^rates|^forwards", parameter) & model == "mascot")) %>%
  mutate(direction = gsub("_[A-Z][A-Za-z_]*$", "", parameter)) %>%
  mutate(direction = ifelse(direction == "rates", "backwards", direction)) %>%
  mutate(parameter = gsub('^[a-z_]*_', '', parameter)) %>%
  left_join(., spatialRates, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol" = "protocol",
                                    "matrix"="matrix", "parameter" ="parameterChar")) %>%
  mutate(value.y = ifelse(is.na(value.y), 0, value.y))

conds <- data[ , c('nSeq', 'protocol', 'matrix')]
conds <- conds[!duplicated(conds),]
regressionDF <- apply(conds, 1, function(x) test_direction(x, migrationRates))
regressionDF <- do.call("rbind", regressionDF)

# Plot
protocols = "bias"
prot = bias
i = 1
ggplot(migrationRates[migrationRates$matrix == i & migrationRates$protocol %in% prot & migrationRates$direction == "forwards", ],
       aes(x = value.y, y = mean, col = direction)) +
  facet_grid(nSeq ~ protocol) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = 'lm') +
  scale_color_manual(values = c(wes_palette("Moonrise2")[2])) +
  theme_minimal() +
  theme(axis.text = element_text(colour = 'black')) +
  # geom_text(data = regressionDF[regressionDF$matrix == i & regressionDF$direction == "backwards" & regressionDF$protocol %in% prot, ],
  #           aes(x = 1e-4, y= 100, label = coef),
  #           col = "black", parse = TRUE,  size = 3, hjust = 0) +
  geom_text(data = regressionDF[regressionDF$matrix == i & regressionDF$direction == "forwards" & regressionDF$protocol %in% prot, ],
            aes(x = 1e-4, y= 40, label = coef),
            col = "black", parse =TRUE, size = 3, hjust = 0) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  labs(x = "Simulated migration events (mean)", y = "Estimated migration events (mean)", col = "Rate direction")

if (protocols == "bias") {
  ggsave(paste0("2.Figures/", runType, "/migration_rates_forwards_", protocols, ".png"), # ,"_m", i, ".png"), 
         width = 13, height = 5)
} else {
  ggsave(paste0("2.Figures/", runType, "/migration_rates_forwards_", protocols, ".png"), # ,"_m", i, ".png"),
         width = 7, height = 5)
}

