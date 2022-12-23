#############################################################
##          ANALYSIS OF LINEAGE MIGRATION COUNTS
#############################################################

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(gridExtra)
library(gtools)
library(ggpubr)

source('R/plot_results.R')

runType = "7demes"
models = c('dta', 'basta', 'mascot', 'glm')
if (runType == "7demes") regions <- regions[1:3]

# Input and output Directory
dirAdjBF = paste0(runType, "/analyses")
dirFigs = paste0("figures/", runType)

# Colors
colModels = c("#00979c", "#883000", "#fa7850", "darkgoldenrod1")
names(colModels) = c("CTMC", "BASTA", "MASCOT", "MASCOT-GLM")

# Plots
allPlots=list()

#############################################################
## Load data
#############################################################
data = data.frame()
allConds = expand.grid(models, c("_150", "_500"))
fileNames = paste0("/lineage_migration_", apply(allConds, 1, paste, collapse = ''),".txt")

for (f in fileNames) {
      fileName = paste0(dirAdjBF, f)
    if (file.exists(fileName)) {
      print(fileName)
      data_tmp = read.table(fileName, stringsAsFactors = F, sep = "\t", header = T)
      if (grepl("_glm_", f)) data_tmp$model = "glm"
      data = bind_rows(data, data_tmp)
    }
}

# Keep only conditions with enough statistical support
selectedRuns = read.table(paste0(dirFigs, "/selected_runs.txt"), 
                          header = T, sep = "\t", stringsAsFactors = F) %>%
  filter(model %in% models) %>%
  select(-matrix)
data = right_join(data, selectedRuns)

# Join with simulated values
simData = read.table(paste0(dirAdjBF, "/lineage_migration_sim_150.txt"), sep = "\t", header = T) %>%
  bind_rows(., read.table(paste0(dirAdjBF, "/lineage_migration_sim_500.txt"), sep = "\t", header = T)) %>%
  select(-model)
data = left_join(data, simData) %>%
  filter(source != destination) %>%
  rename('value.y' = 'value') %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         model = factor(model, levels = models, labels = names(colModels))) 


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
  group_by(nSeq, protocol, model) %>%
  summarise(perc = 100 * sum(value.y >= X2_5_hpd & value.y <= X97_5_hpd) / n() ) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences"))) %>%
  data.frame() 

# Data
data = mutate(data, 
              nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences")))


#############################################################
## Get weighted interval score
#############################################################
# Compute WIS
irs = colnames(data)[grepl("^X[0-9_]+$", colnames(data))]
irs = as.numeric(gsub("_", ".", gsub("X", "", irs)))
irs = irs[irs>50]

wis_data = lapply(irs, function(x) {
  
  X1 = as.character(x)
  X2 = as.character(100-x)
  colNames = c(paste0("X", gsub("\\.", "_", X1)), paste0("X", gsub("\\.", "_", X2)))
  
  data %>%
    mutate(value.y = ifelse(is.na(value.y), 0, value.y)) %>%
    select(all_of(colNames), X50, value.y, source, destination, nSim, nSeq, protocol, model) %>%
    dplyr::rename(point = value.y, smooth_value = X50, upper = colNames[1], lower = colNames[2]) %>%
    mutate(ir = x)
})

#############################################################
## Regression coefficients
#############################################################
n_r = ifelse(grepl("7", runType), 7, 3)
data %>%
  group_by(nSeq) %>%
  summarise(n_mis = sum(value.y == 0), n_tot = n())

# Bias, calibration, mean relative 95% HPD width, average 
# relative error
allPlots = list()
i = 1 
s = c("Mean relative\nbias", "Kendall's tau\n(Correlation)", "Calibration", 
      "Mean relative\n95% HPD width", "Weighted Interval score")


for (statistic in s) {
  
  if (statistic == "Kendall's tau\n(Correlation)") statistic_df = rename(regressionDF, stat = tau)
  if (statistic == "Calibration") statistic_df = rename(calibration, stat = perc)
  if (statistic == "Mean relative\n95% HPD width") statistic_df = data %>%
      filter(value.y > 0) %>%
      group_by(protocol, nSeq, model) %>%
      summarise(stat = mean( (X97_5_hpd - X2_5_hpd) / value.y, na.rm = T) * 100)
  if (statistic == "Mean relative\nbias") statistic_df = data %>%
      filter(value.y > 0) %>%
      group_by(protocol, nSeq, model) %>%
      summarise(stat = mean( (X50 - value.y) / value.y, na.rm = T) * 100 )
  
  if (statistic == "Weighted Interval score") {
    statistic_df = data %>%
      mutate(value.y = ifelse(is.na(value.y), 0, value.y)) %>%
      mutate(stat = WIS(wis_data)) 
    
    # Systematic bias
    allPlots[[i]] = statistic_df %>%
      filter(protocol %in% bias) %>%
      mutate(protocol = factor(protocol, 
                               levels = bias,
                               labels = gsub("biased_", "", bias))) %>%
      ggplot(., aes(x = protocol, y = stat, col = model, group = interaction(model, protocol))) +
      geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
      geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),  
                 alpha = 0.4, size = 0.8, aes(col = as.character(model))) +
      facet_grid(nSeq~.) +
      scale_color_manual(values = colModels) +
      theme_light() +
      theme(panel.background = element_blank()) +
      labs(x = "", y = "WIS", col = "", linetype = "") +
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
                                 labels = c("uni.\nsurv.", "region", "region+\nyear"))) %>%
        filter(surveillance_bias == b) %>%
        ggplot(., aes(x = protocol, y = stat, col = model, group = interaction(model, protocol))) +
        geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
        geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
                   size = 0.8, alpha = 0.4, aes(col = as.character(model))) +
        facet_grid(nSeq~.) +
        scale_color_manual(values = colModels) +
        theme_light() +
        theme(plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
              axis.text.y = element_text(margin = margin(r = 0)),
              panel.background = element_blank()) +
        labs(x = "", y = "") +
        rremove("xlab") +
        rremove("ylab")
      
      k = k+1
    }
  } else {
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
      theme_light()+
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
                                 labels = c("uni.\nsurv.", "region", "region+\nyear"))) %>%
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
  }

  if (statistic == "Kendall's tau\n(Correlation)") {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(c(0,1))
  }
  if (statistic == "Calibration") {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(c(0,100))
  }
  if (statistic %in% c("Mean relative\nbias", "Mean relative\n95% HPD width")) {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] +
        ylim(min(statistic_df$stat),max(statistic_df$stat))
  }
  
  if (!statistic %in% c("Weighted Interval score", "Mean relative\n95% HPD width")) {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + rremove("x.text") + rremove("x.ticks")
  }
  
  i = i+3
}

# All plots
p1 = ggarrange(plotlist = allPlots[c(1:12)], ncol = 3, nrow = 4, 
               labels = c("A", "E", "", "B", "F", "", "C", "G", "", "D", "H", ""), 
               widths = rep(c(1.3,1,1), 4), legend = "top", common.legend = T,
               align = "hv")
  
p2 = ggarrange(plotlist = allPlots[c(13:15)], ncol = 3,
               labels = c("E", "J", ""), widths = c(1.5,1,1), 
               vjust = 1, hjust = 0, align = "hv", legend = "none")

tosave = ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(3.3,1), align = "hv", 
                   common.legend = T, legend = "bottom")

ggsave(paste0("2.Figures/", runType, "/lineage_migration_counts.png"), tosave, width = 9, height = 12)
