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
library(grid)
library(gtable)
library(gridExtra)

source('../R/plot_results.R')

runType = "radiation_bis"
category = paste0('HKY_', runType)

#############################################################
## Load data
#############################################################
## Simulation data
m = data.frame()
for (i in category) {
  fileName = paste0(i, '/analyses/migration_events_sim1to10.txt') 
  if (file.exists(fileName)) {
    mTemp = read.table(fileName, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
    m = bind_rows(m, mTemp) 
  }
}


m = filter(m, protocol != 'stratified') %>%
  #mutate(protocol = gsub("_[0-9]*$", "", protocol)) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance))) %>%
  mutate(parameter = recode(parameter,
                            epidemicDuration = "treeHeight"))

## BEAST data
data = read.table(paste0("2.Figures/", runType, "/2.selected_data.txt"),
                  stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)),
         model = gsub("_.*$", "", model))

## Selected runs
selected = read.table(paste0("2.Figures/", runType, "/2.selected_runs.txt"),
                      stringsAsFactors = FALSE,
                      header = TRUE, sep = "\t")

#############################################################
##  GENETIC PARAMETERS
#############################################################
## Clock rate
clockRateD = data %>%
  filter(parameter == "clockRate")

ggplot(clockRateD, aes(x = protocol, y = X50, fill = model)) +
  facet_grid(nSeq ~ .) +
  geom_boxplot() +
  scale_fill_manual(values = c(wes_palette("Darjeeling1")[c(2,3)])) +
  theme_minimal() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = clockRate) +
  labs(x = "", y = "Clock rate (/site/year)", col = "Sampling protocol") 

## Kappa and frequencies
site = data %>%
  filter(parameter %in% c('freqA', 'freqC', 'freqG', 'freqT', 'kappa'))

siteParam = data.frame(X50 = c(kappa, freqA, freqC, freqG, freqT), 
                       parameter =  c("kappa", "freqA", "freqC", "freqG", "freqT"))

ggplot(site, aes(x = protocol, y = X50, fill = model)) +
  facet_grid(parameter ~ ., scales = "free_y") +
  geom_boxplot() +
  scale_fill_manual(values = c(wes_palette("Darjeeling1")[c(2,3)])) +
  theme_minimal() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45,  hjust = 1)) + 
  geom_hline(data = siteParam, aes(yintercept = X50)) +
  labs(x = "", y = "", col = "Sampling protocol")

# Root age 
rootAgeS <- m %>%
  filter(parameter == "treeHeight") %>%
  mutate(value = as.numeric(value))

rootAge <-  data %>%
  filter(parameter == "treeHeight") %>%
  left_join(., rootAgeS, by = c("parameter" = "parameter", "nSim" = "nSim", 
                                "nSeq" = "nSeq", "matrix"="matrix", "protocol"="protocol"))

conds <- data[, c('nSeq', 'protocol')]
conds <- conds[!duplicated(conds),]
regressionAge = apply(conds, 1, function(x) test_eqn(x, rootAge))
regressionAge = do.call("rbind", regressionAge)

ggplot(rootAge[rootAge$protocol %in% bias, ], 
       aes(x = value.y, y = X50, ymin = X2_5, ymax = X97_5, col = model)) +
  facet_grid(nSeq ~ protocol) +
  geom_pointrange(fatten = 1.5) +
  scale_color_manual(values = c(wes_palette("Darjeeling1")[2:3])) +
  geom_abline(slope = 1, intercept = 0, col = "grey50", linetype = "longdash") +
  geom_text(data = regressionAge[regressionAge$protocol %in% bias & regressionAge$model == "dta", ],
            aes(x = 19, y=30, label = coef),
            hjust = 0, col = "black", parse = TRUE, size = 3, inherit.aes = FALSE) +
  # geom_text(data = regressionAge[regressionAge$protocol %in% bias & regressionAge$model == "mascot", ],
  #           aes(x = 19, y=30, label = coef),
  #           hjust = 0, col = "black", parse = TRUE, size = 3, inherit.aes = FALSE) +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(colour = 'black')) + 
  labs(x = "Simulated tree height (years)", y = "Inferred tree height (years, median and 95% CI)", col = "BEAST model") 
ggsave(paste0("2.Figures/", runType, "/rootage_bias.png"), width = 13, height = 5)

ggplot(rootAge[rootAge$protocol %in% surveillance, ], 
       aes(x = value.y, y = X50, ymin = X2_5, ymax = X97_5, col = model)) +
  facet_grid(nSeq ~ protocol) +
  geom_pointrange(fatten = 1.5) +
  scale_color_manual(values = c(wes_palette("Darjeeling1")[2:3])) +
  geom_abline(slope = 1, intercept = 0, col = "grey50", linetype = "longdash") +
  geom_text(data = regressionAge[regressionAge$protocol %in% surveillance & regressionAge$model == "dta", ],
            aes(x = 19, y=30, label = coef),
            hjust = 0, col = "black", parse = TRUE, size = 3, inherit.aes = FALSE) +
  # geom_text(data = regressionAge[regressionAge$protocol %in% surveillance & regressionAge$model == "mascot", ],
  #           aes(x = 19, y=30, label = coef),
  #           hjust = 0, col = "black", parse = TRUE, size = 3, inherit.aes = FALSE) +
  theme_minimal() +
  theme(axis.text = element_text(colour = 'black'), legend.position = "none") + 
  labs(x = "Simulated tree height (years)", y = "Inferred tree height (years, median and 95% CI)", col = "BEAST model") 
ggsave(paste0("2.Figures/", runType, "/rootage_surveillance.png"), width = 6, height = 5)

# RMSE Tables
rmseTab = 
  data %>%
  filter(parameter %in% c('freqA', 'freqC', 'freqG', 'freqT', 'kappa', 'clockRate', 'treeHeight')) %>%
  left_join(., rootAgeS, by = c("parameter" = "parameter", "nSim" = "nSim","nSeq" = "nSeq", 
                                "matrix"="matrix", "protocol"="protocol")) %>%
  group_by(nSeq, protocol, matrix, model, parameter) %>%
  summarise(rmse = RMSE(X50, value.y, parameter)) %>%
  mutate(rmse = ifelse(parameter %in% c('freqA', 'freqC', 'freqG', 'freqT'), round(rmse, 5), round(rmse, 3))) %>%
  spread(., model, rmse) %>%
  data.frame()

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.9)),
  colhead = list(fg_params=list(cex = 1)),
  rowhead = list(fg_params=list(cex = 1)))

for (i in c('freqA', 'freqC', 'freqG', 'freqT', 'kappa', 'clockRate', 'treeHeight')) {
  
  tab <- rmseTab[rmseTab$parameter == i, ]
  tab <- tab[order(tab$protocol, tab$nSeq), colnames(tab) != "parameter"]
  
  
  g = tableGrob(tab, rows = NULL, theme = mytheme)
  g = gtable_add_grob(g,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 2, b = nrow(g), l = 1, r = ncol(g))
  g = gtable_add_grob(g,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, b = nrow(g), l = 1, r = ncol(g))
  
  png(paste0("2.Figures/", runType, "/rmse_", i,".png"))
  grid.draw(g)
  dev.off()
}  

# Verify that ESS of genetic parameters are greater than
# threshold in all included runs
data %>% 
  select(model, nSeq, nSim, protocol, matrix, parameter, ESS) %>%
  filter(parameter %in% c('freqA', 'freqC', 'freqG', 'freqT', 'kappa', 'clockRate', 'treeHeight')) %>%
  filter(ESS < 200) %>%
  spread(parameter, ESS)

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
  filter((grepl("^bssvs_", parameter) & model == "mascot")) %>% 
  mutate(direction = gsub("^bssvs_|_[A-Z].*$", "", parameter)) %>%
  # filter(BF >= 3) %>%
  # filter((grepl("^rates|^forwards", parameter) & model == "mascot")) %>%
  # mutate(direction = gsub("_.*$", "", parameter)) %>%
  mutate(direction = ifelse(direction == "rates", "backwards", direction)) %>%
  mutate(parameter = gsub('^[a-z_]*_', '', parameter)) %>%
  left_join(., spatialRates, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol" = "protocol",
                                    "matrix"="matrix", "parameter" ="parameterChar")) %>%
  mutate(value.y = ifelse(is.na(value.y), 0, value.y))

conds <- data[ , c('nSeq', 'protocol', 'matrix')]
conds <- conds[!duplicated(conds),]
regressionDF <- apply(conds, 1, function(x) test_direction(x, migrationRates))
regressionDF <- do.call("rbind", regressionDF)
regressionDF$protocol = factor(regressionDF$protocol, levels = c(bias, surveillance))

# Plot
protocols = "surveillance"
prot = surveillance
i = 1

ggplot(migrationRates[migrationRates$protocol %in% prot, ],
       aes(x = value.y, y = mean, col = direction)) +
  facet_grid(nSeq ~ protocol) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = 'lm') +
  scale_color_manual(values = c(wes_palette("Moonrise2")[c(1,2)])) +
  theme_minimal() +
  theme(axis.text = element_text(colour = 'black')) +
  geom_text(data = regressionDF[regressionDF$direction == "backwards" & regressionDF$protocol %in% prot, ],
            aes(x = 1e-4, y= 4, label = coef),
            col = "black", parse = TRUE,  size = 3, hjust = 0) +
  geom_text(data = regressionDF[regressionDF$direction == "forwards" & regressionDF$protocol %in% prot, ],
            aes(x = 1e-4, y= 10, label = coef),
            col = "black", parse =TRUE, size = 3, hjust = 0) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  labs(x = "Simulated migration events (mean)", y = "Estimated migration events (mean)", col = "Rate direction")

if (protocols == "bias") {
  ggsave(paste0("2.Figures/", runType, "/migration_rates_", protocols, ".png"), # ,"_m", i, ".png"), 
         width = 13, height = 5)
} else {
  ggsave(paste0("2.Figures/", runType, "/migration_rates_", protocols, ".png"), # ,"_m", i, ".png"),
         width = 7, height = 5)
}


#############################################################
##  DTA model
#############################################################
# Root location - cropped tree 
rootLocS <- m %>%
  filter(parameter == "mrcaLocation") %>%
  mutate(parameter = gsub(" ", "", value)) %>%
  dplyr::select(-value, -parameterChar)

rootLoc <- 
  data %>%
  filter(grepl("root_", parameter)) %>%
  mutate(parameter = gsub('root_', '', parameter)) %>%
  right_join(., rootLocS, by = c("nSim", "nSeq","matrix", "protocol", "parameter"))

ggplot(rootLoc[rootLoc$protocol %in% bias, ], aes(x = protocol, y = value, fill = model)) +
  facet_grid(nSeq ~ ., labeller = label_both) +
  geom_boxplot(width = 0.3) +
  scale_fill_manual(values = c(wes_palette("Darjeeling1")[2])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Sampling protocol", y = "Root state posterior probability", fill = "BEAST model")
ggsave(paste0("2.Figures/", runType, "/rootlocation_cropped_tree_bias.png"), width = 5, height = 5)

## Migration rates 
# conds <- data[ , c('nSim', 'nSeq', 'protocol', 'matrix')] 
# conds <- conds[!duplicated(conds),]
# 
# spatialRates <- m %>%
#   filter(parameter %in% c("locationFrequency", "migrationRate")) %>%
#   mutate(value = as.numeric(value)) %>%
#   mutate(source = gsub("_[a-zA-Z]*$", "", parameterChar),
#          dest = gsub("^[a-zA-Z]*_", "", parameterChar))
# 
# spatialRatesNormalized = apply(conds, 1, function(x) normalizeRates(x, spatialRates))
# spatialRatesNormalized = do.call("rbind", spatialRatesNormalized)
# 
# migrationRates <- data %>%
#   filter(grepl("bssvs_", parameter)) %>%
#   filter(BF >= 3) %>%
#   mutate(parameter = gsub('^[a-z_]*_', '', parameter)) %>%
#   left_join(., spatialRatesNormalized, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol" = "protocol",
#                                               "matrix"="matrix", "parameter" ="parameterChar")) %>%
#   mutate(value.y = ifelse(is.na(value.y), 0, value.y))
# 
# conds <- data[ , c('nSeq', 'protocol', 'matrix')] 
# conds <- conds[!duplicated(conds),]
# regressionDF <- apply(conds, 1, function(x) test_eqn(x, migrationRates, stat = "mean"))
# regressionDF <- do.call("rbind", regressionDF)

spatialRates <- m %>%
  filter(parameter %in% c("locationFrequency", "migrationRate")) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(source = gsub("_[a-zA-Z]*$", "", parameterChar),
         dest = gsub("^[a-zA-Z]*_", "", parameterChar))

migrationRates <- data %>%
  filter((grepl("^bssvs_", parameter) & model == "dta"), BF >= 3) %>%
  mutate(parameter = gsub('^[a-z_]*_', '', parameter)) %>%
  left_join(., spatialRates, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol" = "protocol",
                                    "matrix"="matrix", "parameter" ="parameterChar")) %>%
  mutate(value.y = ifelse(is.na(value.y), 0, value.y))

conds <- data[ , c('nSeq', 'protocol', 'matrix')]
conds <- conds[!duplicated(conds),]
regressionDF <- apply(conds, 1, function(x) test_eqn(x, migrationRates, stat = "mean"))
regressionDF <- do.call("rbind", regressionDF)
regressionDF$protocol <- factor(regressionDF$protocol, levels = c(bias, surveillance))

# Plot
protocols = "bias"
prot = bias

ggplot(migrationRates[migrationRates$protocol %in% prot, ], 
       aes(x = value.y, y = mean, col = model)) +
  facet_grid(nSeq ~ protocol) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = 'lm') +
  scale_color_manual(values = c(wes_palette("Darjeeling1")[2])) +
  theme_minimal() + 
  theme(axis.text = element_text(colour = 'black'), legend.position = "none") +
  geom_text(data = regressionDF[regressionDF$model == "dta" & regressionDF$protocol %in% prot, ],
            aes(x = 1e-4, y=4, label = coef),
            col = "black",  parse  = TRUE, size = 3, hjust = 0) +
  # scale_x_continuous(trans='log10') +
  # scale_y_continuous(trans='log10') +
  labs(x = "Simulated migration events (mean)", y = "Estimated migration events (mean)", col = "BEAST model")

if (protocols == "bias") {
  ggsave(paste0("2.Figures/", runType, "/migration_rates_bssvs_", protocols, ".png"), # ,"_m", i, ".png"), 
         width = 13, height = 5)
} else {
  ggsave(paste0("2.Figures/", runType, "/migration_rates_bssvs_", protocols, ".png"), # ,"_m", i, ".png"),
         width = 7, height = 5)
}

# Markov jumps
spatialCounts <- m %>%
  filter(parameter == "migrationCount") %>%
  mutate(value = as.numeric(value))

migrationEvents <- data %>%
  filter(grepl("nMigration_", parameter)) %>%
  mutate(parameter = gsub('nMigration_', '', parameter)) %>%
  left_join(., spatialCounts, by = c("nSim" = "nSim", "nSeq" = "nSeq", 
                                     "matrix"="matrix", "protocol"="protocol", "parameter" ="parameterChar"))

for (n in 1:2) {
  ggplot(migrationEvents[migrationEvents$protocol %in% bias &
                           !is.na(migrationEvents$value.y) & !is.na(migrationEvents$X50), ], 
         aes(x = value.y, y = X50, col = model)) +
    facet_grid(nSeq ~ protocol, labeller = label_both) +
    geom_point() +
    geom_smooth(formula = 'y ~ x', method='lm') +
    geom_abline(slope = 1, intercept = 0, col = "grey50", linetype = "longdash") +
    scale_color_manual(values = c(wes_palette("Darjeeling1")[2])) +
    theme_minimal() + 
    theme(legend.position = "none", axis.text = element_text(colour = "black")) +
    labs(x = "Simulated migration events", y = "Estimated migration events")
  ggsave(paste0("2.Figures/migration_events_bias_m", n, ".png"), width = 11, height = 5)
} 

for (n in 1:2) {
  ggplot(migrationEvents[migrationEvents$protocol %in% surveillance &
                           !is.na(migrationEvents$X50) & !is.na(migrationEvents$value.y), ], 
         aes(x = value.y, y = X50, col = model)) +
    facet_grid(nSeq ~ protocol) +
    geom_point() +
    geom_smooth(formula = 'y ~ x', method='lm') +
    geom_abline(slope = 1, intercept = 0, col = "grey50", linetype = "longdash") +
    scale_color_manual(values = c(wes_palette("Darjeeling1")[2])) +
    theme_minimal() + 
    theme(axis.text = element_text(color = "black"), legend.position = "none") +
    labs(x = "Simulated migration events", y = "Estimated migration events")
  ggsave(paste0("2.Figures/migration_events_surveillance_m", n, ".png"), width = 6, height = 5)
} 
