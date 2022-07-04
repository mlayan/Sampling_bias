#############################################################
##                ANALYSIS OF INTRODUCTION DATES
#############################################################

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(gridExtra)
library(gtools)
library(ggpubr)

source('R/plot_results.R')

runType = "3demes"
category = paste0('HKY_', runType)
models = c("mascot", 'dta', 'basta', "glm")
if (runType == "3demes") regions <- regions[1:3]

# Input and output Directory
dirFigs = ifelse(exists("subcategory"),
                 paste0("2.Figures/", runType, "/", subcategory),
                 paste0("2.Figures/", runType))
dirAdjBF = paste0(category, "/analyses")

# Colors
colModels = c("#00979c", "#883000", "#fa7850", "darkgoldenrod1")
names(colModels) = c("CTMC", "BASTA", "MASCOT", "MASCOT-GLM")

# Plots
allPlots = list()

#############################################################
## Load data
#############################################################
data = data.frame()

allConds = expand.grid(models, c(150,500))
allConds = apply(allConds, 1, paste, collapse = "_")
files = paste0("/introduction_dates_", allConds, ".txt")

for (f in files) {
  if (!grepl("dta_500", f)) {
      fileName = paste0(dirAdjBF, f)
    if (file.exists(fileName)) {
      print(fileName)
      if (nrow(data) == 0) {
        data = read.table(fileName, stringsAsFactors = F, sep = "\t", header = T)
      } else {
        data = bind_rows(data, 
                         read.table(fileName, stringsAsFactors = F, sep = "\t", header = T))
      } 
    }
  }
}

# Keep only conditions with enough statistical support
selectedRuns = read.table(paste0(dirFigs, "/2.selected_runs.txt"), 
                          header = T, sep = "\t", stringsAsFactors = F) %>%
  filter(model %in% models) %>%
  select(-matrix)
data = right_join(data, selectedRuns)

# Load simulated parameters
simData1 = read.table(paste0(dirAdjBF, "/introduction_dates_sim_150.txt"), sep = "\t", header = T)
simData = rbind(
  simData1, 
  read.table(paste0(dirAdjBF, "/introduction_dates_sim_500.txt"), sep = "\t", header = T)
) %>%
  select(-model)

# Join estimated and simulated parameters
data = left_join(data, simData) %>%
  rename('value.y' = 'value') %>%
  filter(!is.na(X50)) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         model = factor(model, levels = c("dta", "basta", "mascot", "glm"), 
                        labels = c("CTMC", "BASTA", "MASCOT", "MASCOT-GLM"))) 

# Number of parameters
nParameters = data %>%
  group_by(nSeq, protocol, model) %>%
  summarise(nrates = n())

#############################################################
## Get correlation and calibration
#############################################################
# Correlation
conds <- data[ , c('nSeq', 'protocol')] 
conds <- conds[!duplicated(conds),]
regressionDF <- apply(conds, 1, function(x) test_eqn(x, data, test = "kendall"))
regressionDF <- do.call("rbind", regressionDF) %>%
  filter(!is.na(p)) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences"))) 

# Calibration
calibration = data %>%
  filter(value.y >= X2_5_hpd, value.y <= X97_5_hpd) %>%
  group_by(nSeq, protocol, model) %>%
  summarise(withinInterval = n()) %>%
  full_join(., nParameters) %>%
  mutate(perc = withinInterval/nrates*100) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences"))) %>%
  data.frame() 

# Data
data = mutate(data, 
              nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences")))

#############################################################
## Regression coefficients
#############################################################
# Bias, calibration, mean relative 95% HPD width, average 
# relative error
allPlots = list()
i = 1 
s = c("Kendall's tau", "Calibration", "Mean relative 95% HPD\nwidth", "Mean relative bias")


for (statistic in s) {
  
  if (statistic == "Kendall's tau") statistic_df = rename(regressionDF, stat = tau)
  if (statistic == "Calibration") statistic_df = rename(calibration, stat = perc)
  if (statistic == "Mean relative 95% HPD\nwidth") statistic_df = data %>%
      group_by(protocol, nSeq, model) %>%
      summarise(stat = mean( (X97_5_hpd - X2_5_hpd) / value.y, na.rm = T))
  if (statistic == "Mean relative bias") statistic_df = data %>%
      group_by(protocol, nSeq, model) %>%
      summarise(stat = mean( (X50 - value.y) / value.y, na.rm = T) )
  
  # Systematic bias
  allPlots[[i]] = statistic_df %>%
    filter(protocol %in% bias) %>%
    mutate(protocol = factor(protocol, 
                             levels = bias,
                             labels = gsub("biased_", "", bias))) %>%
    ggplot(., aes(x = protocol, y = stat, col = model, shape = nSeq, group = interaction(model, nSeq))) +
    geom_point(size = 2) +
    geom_line(size = 0.8) +
    scale_color_manual(values = colModels) +
    scale_shape_manual(values = c(16,17)) +    
    theme_light() +
    theme(panel.background = element_blank()) +
    labs(x = "", y = statistic, col = "", shape = "") +
    rremove("xlab")
  
  # Surveillance bias
  k = 1
  for (b in c(10,20)) {
    
    allPlots[[i+k]] = statistic_df %>% 
      filter(protocol %in% surveillance) %>%
      separate(protocol, 
               into = c("protocol", "surveillance_bias"), 
               sep = "_") %>%
      mutate(protocol = factor(protocol, 
                               levels = c("uniformS", "maxPerRegion", "maxPerRegionYear"), 
                               labels = c("uni. surv.", "region", "region+year"))) %>%
      filter(surveillance_bias == b) %>%
      ggplot(., aes(x = protocol, y = stat, col = model, shape = nSeq, group = interaction(nSeq, model))) +
      geom_point(size = 2) +
      geom_line(size = 0.8) +
      scale_color_manual(values = colModels) +
      scale_shape_manual(values = c(16,17)) +
      theme_light() +
      theme(plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
            axis.text.y = element_text(margin = margin(r = 0)),
            panel.background = element_blank()) +
      labs(x = "", y = "") +
      rremove("xlab") + 
      rremove("ylab")
    
    k = k+1
  }
  
  if (statistic == "Kendall's tau") {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(-0.5,1)
  }
  if (statistic == "Calibration") {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(0,100)
  }
  if (statistic %in% c("Mean relative 95% HPD width", "Mean relative bias")) {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(min(statistic_df$stat), max(statistic_df$stat))
  }
  
  if (i < (length(s)-1)*3) {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + rremove("x.text") + rremove("x.ticks")
  }
  
  i = i+3
}

arrange = ggarrange(plotlist = allPlots, ncol = 3, nrow = 4,
                    labels = c("A", "E", "", "B", "F", "", "C", "G", "", "D", "H", ""), 
                    widths = rep(c(1.3,1,1), 4), vjust = 1, hjust = 0,
                    align = "hv", common.legend = T, legend = "bottom")
ggsave(paste0("2.Figures/", runType, "/introduction_dates.png"), arrange, width = 9, height = 9)

