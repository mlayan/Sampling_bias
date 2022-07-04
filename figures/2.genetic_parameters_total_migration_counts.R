#############################################################
##        COMPARISON OF DTA AND APPROXIMATIONS OF 
##                  STRUCTURED COALESCENT
##
## Maylis Layan
#############################################################

rm(list = ls())
library(tidyverse)
library(wesanderson)
library(grid)
library(gtable)
library(gridExtra)
library(gtools)
library(ggpubr)

source('R/plot_results.R')

runType = "3demes"
ftl = list(size = 20, color = "black", face = "bold")

category = paste0('HKY_', runType)
if (runType == "3demes") regions = regions[1:3]
fluxes <- apply(permutations(length(regions), 2, regions), 1, function(x) paste0(x[1], "_", x[2]))

# Input and output Directory
directory = paste0("figures/", runType)

# Colors
colModels = c("#00979c", "#883000", "#fa7850", "darkgoldenrod1")
names(colModels) = c("CTMC", "BASTA", "MASCOT", "MASCOT-GLM")


#############################################################
## Load data
#############################################################
## Simulation data
m = data.frame()
for (i in category) {
  fileName = paste0(i, '/analyses/migration_events.txt')
  
  if (file.exists(fileName)) {
    mTemp = read.table(fileName, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
    m = bind_rows(m, mTemp) 
  }
}

m = m %>% 
  filter(protocol %in% c(bias, surveillance)) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)))

if (runType == "7demes") m$matrix = runType

## BEAST data
data = read.table(paste0(directory, "/2.selected_data.txt"),
                  stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>%
  filter(!(model == "basta" & nSeq == 500)) %>% 
  mutate(model = factor(model, levels = c("dta", "basta", "mascot", "glm"), labels = names(colModels))) 

## Selected runs
selected = read.table(paste0(directory, "/2.selected_runs.txt"),
                      stringsAsFactors = FALSE,
                      header = TRUE, sep = "\t")


## Tree topology
topologies = read.table(paste0(category, '/analyses/topology_regression.txt'), 
                        stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>%
  right_join(., selected, by = c("model", "nSim", "protocol", "nSeq", "matrix")) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         model = factor(model, 
                        c("dta", "basta", "mascot", "glm"), 
                        c("CTMC", "BASTA",  "MASCOT", "MASCOT-GLM"))) %>%
  mutate(nSeq = factor(nSeq, c(150, 500), c("150 sequences", "500 sequences")))

#############################################################
##  GENETIC PARAMETERS
#############################################################
## Figure 2--------------------------------------------------
figure2 = list()

figure2[[1]] = topologies %>%
  filter(protocol %in% bias) %>%
  mutate(protocol = factor(protocol, bias, c("uniform", "2.5", "5", "10", "20", "50"))) %>%
  ggplot(.,aes(x = protocol, y = coef, group = interaction(protocol, model))) +
  facet_grid(nSeq ~ .) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),  
             alpha = 0.4, aes(col = as.character(model))) +
  scale_color_manual(values = colModels) +
  theme_light() +
  labs(x = "", y = "R2 of linear regression on tree topology", col = "") 

clockRateD = data %>%
  filter(parameter == "clockRate") %>%
  mutate(nSeq = factor(nSeq, levels = c(150, 500), labels = c("150 sequences", "500 sequences")), 
         protocol = factor(protocol, levels = c(bias, surveillance)))

figure2[[2]] = clockRateD %>%
  filter(protocol %in% bias) %>%
  mutate(protocol = factor(protocol, bias, c("uniform", "2.5", "5", "10", "20", "50"))) %>%
  ggplot(.,aes(x = protocol, y = X50, group = interaction(protocol, model))) +
  facet_grid(nSeq ~ .) +
  geom_hline(yintercept = clockRate, size = 1) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),  
             alpha = 0.4, aes(col = as.character(model))) +
  scale_color_manual(values = colModels) +
  theme_light() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(x = "", y = "Evolutionary rate (median, /site/year)", col = "") 

arranged = ggarrange(plotlist = figure2, ncol = 1, nrow = 2, 
                     labels = c("B", "C"), common.legend = T, legend = "bottom", align = "hv")
ggsave(paste0("2.Figures/", runType, "/figure2.png"), arranged, width = 7, height = 8)  
 

## Tree topology--------------------------------------------------------
treeTop = list()
ymin = min(topologies$coef, na.rm = T) 

treeTop[[1]] = topologies %>%
  filter(protocol %in% bias) %>%
  mutate(protocol = factor(protocol, bias, c("uniform", "2.5", "5", "10", "20", "50"))) %>%
  ggplot(.,aes(x = protocol, y = coef, group = interaction(protocol, model))) +
  facet_grid(rows = vars(nSeq)) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),  
             alpha = 0.4, aes(col = as.character(model))) +
  scale_color_manual(values = colModels) +
  theme_light() +
  theme(legend.position = "none") +
  ylim(ymin, 1) +
  labs(x = "", y = "R2 of linear regression on tree topology", col = "") +
  rremove("xlab")
  

k=1
for (b in c(10,20)) {
  treeTop[[k+1]] = topologies %>%
    filter(protocol %in% surveillance) %>%
    separate(protocol, into = c("protocol", "surveillance_bias"), sep = "_") %>%
    filter(surveillance_bias == b) %>%
    mutate(protocol = factor(protocol, 
                             levels = c("uniformS", "maxPerRegion", "maxPerRegionYear"), 
                             labels = c("uni. surv.", "region", "region+year"))) %>%
    ggplot(.,aes(x = protocol, y = coef, group = interaction(protocol, model))) +
    facet_grid(rows = vars(nSeq)) +
    geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),  
               alpha = 0.4, aes(col = as.character(model))) +
    scale_color_manual(values = colModels) +
    theme_light() +
    labs(x = "", y = "R2 of linear regression on tree topology", col = "") +
    ylim(ymin, 1)+
    rremove("xlab")
  
  if (b == 10) treeTop[[1+k]] = treeTop[[1+k]] + labs(x = "", y = "R2 of linear regression on tree topology", col = "")
  if (b == 20) treeTop[[1+k]] = treeTop[[1+k]] + labs(x = "", y = "", col = "")
  
  k = k+1
}


arranged = ggarrange(treeTop[[1]], 
                     ggarrange(plotlist = treeTop[2:3], 
                               common.legend = T, legend = "bottom", ncol = 2, hjust = 0, align = "v", labels = c("B", "C")), 
                     nrow = 2, hjust = 0, labels = c("A", ""))

ggsave(paste0("2.Figures/", runType, "/topologies.png"), arranged, width = 8, height = 8)
  

## Kappa and frequencies------------------------------------------------
p_name = c('Base frequency A', 'Base frequency C', 'Base frequency G', 
           'Base frequency T', 'Transition-transversion ratio', 'Evolutionary rate')
names(p_name) = c('freqA', 'freqC', 'freqG', 'freqT', 'kappa', 'clockRate')

for (p in names(p_name)) {
  geneticParams = list()
  
  ymin = min(data$X50[data$parameter == p]) 
  ymax = max(data$X50[data$parameter == p])
  
  geneticParams[[1]] = data %>%
    filter(parameter %in% p, protocol %in% bias) %>%
    mutate(
      protocol = factor(protocol, bias, c("uniform", "2.5", "5", "10", "20", "50")),
      nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences"))) %>%
    ggplot(., aes(x = protocol, y = X50, group=interaction(protocol, model))) +
    facet_grid(rows = vars(nSeq)) +
    geom_hline(yintercept = get(p)) +
    geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.4, aes(col = as.character(model))) +
    scale_color_manual(values = colModels) +
    theme_light() +
    theme(legend.position = "none") +
    ylim(c(ymin, ymax)) +
    labs(x = "", y = p_name[p]) +
    rremove("xlab")
  
  k=1
  for (b in c(10,20)) {
    geneticParams[[1+k]] = data %>%
      filter(parameter == p, protocol %in% surveillance) %>%
      separate(protocol, into = c("protocol", "surveillance_bias"), sep = "_") %>%
      filter(surveillance_bias == b) %>%
      mutate(protocol = factor(protocol, 
                               levels = c("uniformS", "maxPerRegion", "maxPerRegionYear"), 
                               labels = c("uni. surv.", "region", "region+year")),
             nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences"))) %>%
      ggplot(., aes(x = protocol, y = X50, group=interaction(protocol, model))) +
      facet_grid(rows = vars(nSeq)) +
      geom_hline(yintercept = get(p)) +
      geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
      geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
                 alpha = 0.4, aes(col = as.character(model))) +
      scale_color_manual(values = colModels) +
      theme_light() +
      ylim(c(ymin, ymax)) +
      rremove("xlab")
    
    if (b == 10) geneticParams[[1+k]] = geneticParams[[1+k]] + labs(x = "", y = p_name[p], col = "")
    if (b == 20) geneticParams[[1+k]] = geneticParams[[1+k]] + labs(x = "", y = "", col = "")
    
    k = k+1
  }
  
  
  arranged = ggarrange(geneticParams[[1]], 
                       ggarrange(plotlist = geneticParams[2:3], common.legend = T, legend = "bottom", ncol = 2, hjust = 0, align = "v", labels = c("B", "C")), 
                       nrow = 2, hjust = 0, labels = c("A", ""))
  
  ggsave(paste0("2.Figures/", runType, "/", p, ".png"), arranged, width = 8, height = 8)
  
}

#############################################################
## ROOT LOCATION
#############################################################
rootPlot = list()

# Data  
rootLocS <- m %>%
  filter(parameter == "mrcaLocation") %>%
  mutate(parameter = gsub(" ", "", value)) %>%
  dplyr::select(-value, -parameterChar)

rootLoc <- data %>%
  filter(grepl("root_", parameter)) %>%
  mutate(parameter = gsub('root_', '', parameter)) %>%
  right_join(., rootLocS, by = c("nSim", "nSeq","matrix", "protocol", "parameter")) %>%
  filter(!is.na(model)) %>%
  mutate(nSeq = factor(nSeq, c(150, 500), c("150 sequences", "500 sequences")))


rootPlot[[1]] = rootLoc %>%
  filter(protocol %in% bias) %>%
  mutate(protocol = factor(protocol, bias, c("uniform", "2.5", "5", "10", "20", "50"))) %>%
  ggplot(., aes(x = protocol, y = value, group = interaction(protocol, model))) +
  facet_grid(nSeq ~ .) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
             alpha = 0.4, aes(col = as.character(model))) +
  theme_light() +
  scale_color_manual(values = colModels) +
  theme(legend.position = "none") +
  labs(x = "", y = "Posterior probability", col = "") +
  ylim(c(0,1))

# Plot surveillance samples
rootLocSurv = rootLoc %>% 
  filter(protocol %in% surveillance) %>%
  separate(protocol, 
           into = c("protocol", "surveillance_bias"), 
           sep = "_") %>%
  mutate(protocol = factor(protocol, 
                           levels = c("uniformS", "maxPerRegion", "maxPerRegionYear"), 
                           labels = c("uni. surv.", "region", "region+year"))) 

i = 2
for (b in c(10,20)) {
  rootPlot[[i]] = rootLocSurv %>%
    filter(surveillance_bias == b) %>%
    ggplot(., aes(x = protocol, y = value, group = interaction(protocol, model))) +
    facet_grid(nSeq ~ .) +
    geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
               alpha = 0.4, aes(col = as.character(model))) +
    theme_light() +
    scale_color_manual(values = colModels) +
    ylim(c(0,1))
  
  if (b == 10) rootPlot[[i]] = rootPlot[[i]] + labs(x = "", y = "Posterior probability", col = "")
  if (b == 20) rootPlot[[i]] = rootPlot[[i]] + labs(x = "", y = "", col = "")
  i = i+1
}

arranged = ggarrange(rootPlot[[1]], 
  ggarrange(plotlist = rootPlot[2:3], labels = c("B", "C"), 
            hjust = 0,align = "v", common.legend = T, legend="bottom"), 
  labels = "A", nrow = 2, hjust = 0)
ggsave(paste0("2.Figures/", runType, "/rootlocation.png"), arranged, width = 8, height = 8)


#############################################################
## Markov jumps
#############################################################
MJPlots = list()

# Simulated total migration counts---------------------------
spatialCounts <- m %>%
  filter(parameter == "migrationCount") %>%
  mutate(value = as.numeric(value))

# Merge simulated and estimated parameters-----------------------
migrationEvents <- data %>%
  filter((grepl("bssvs_nMigration_", parameter) & BF >=3) | (grepl("nMigration_", parameter) & model == "MASCOT-GLM")) %>%
  mutate(parameter = gsub('.*nMigration_', '', parameter)) %>%
  left_join(., spatialCounts, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol"="protocol", 
                                     "parameter" ="parameterChar")) %>%
  mutate(value.y = ifelse(is.na(value.y), 0, value.y)) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance))) 

conds <- migrationEvents[ , c('nSeq', 'protocol')] 
conds <- conds[!duplicated(conds),]
regressionMJ <- apply(conds, 1, function(x) test_eqn(x, migrationEvents, test = "kendall"))

regressionMJ <- do.call("rbind", regressionMJ) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)), 
         nSeq = factor(nSeq, c(150, 500), c("150 sequences", "500 sequences")))

migrationEvents <- mutate(migrationEvents, 
                          nSeq = factor(nSeq, c(150,500), c("150 sequences", "500 sequences")))

migrationEventsCalibration = migrationEvents %>%
  group_by(model, protocol, nSeq) %>%
  summarise(perc = 100 * sum((BF >= 3 | is.na(BF)) & value.y >= X2_5_hpd & value.y <= X97_5_hpd) / sum(BF >= 3 | is.na(BF)))


# Supplementary Figure 8-----------------------------------
y1 = max(migrationEvents$X97_5_hpd[migrationEvents$model == "CTMC" & migrationEvents$protocol %in% bias])

labels = regressionMJ %>%
  filter(protocol %in% bias, model== "CTMC") %>%
  mutate(protocol = gsub("biased_", "", protocol)) %>%
  mutate(protocol = factor(protocol, c('uniform', '2.5', '5', '10', '20', '50')))

p = migrationEvents %>%
  filter(protocol %in% bias, model == "CTMC") %>%
  mutate(protocol = gsub("biased_", "", protocol)) %>%
  mutate(protocol = factor(protocol, c('uniform', '2.5', '5', '10', '20', '50'))) %>%
  ggplot(., aes(x = value.y, y = X50, ymin = X2_5_hpd, ymax = X97_5_hpd, col = model)) +
  geom_pointrange(fatten = 0.5) +
  geom_abline(slope = 1, intercept = 0, col = "grey20", linetype = "dashed") +
  geom_smooth(formula = 'y ~ x', method='lm') +
  facet_grid(rows = vars(nSeq), cols=vars(protocol)) +
  geom_text(data = labels,aes(x = 0, y=y1, label = coef), inherit.aes = FALSE,
            col = "black", size = 3, hjust = 0, parse=T) +
  scale_color_manual(values = colModels[md]) +
  theme_light() + 
  theme(legend.position = "none") +
  labs(x = "Simulated migration events", 
       y = "Estimated migration events (median, 95%-HPD)")

ggsave(paste0(directory, "/supp_figure8.png"), p, width = 10, height = 5)


# Plot for the paper---------------------------------------
# Bias, calibration, mean relative 95% HPD width, average 
# relative error, WIS

# Compute WIS
irs = colnames(data)[grepl("^X[0-9_]+$", colnames(data))]
irs = as.numeric(gsub("_", ".", gsub("X", "", irs)))
irs = irs[irs>50]

wis_data = lapply(irs, function(x) {
  
  X1 = as.character(x)
  X2 = as.character(100-x)
  colNames = c(paste0("X", gsub("\\.", "_", X1)), paste0("X", gsub("\\.", "_", X2)))
  
  data %>%
    filter((grepl("bssvs_nMigration_", parameter) & BF >=3) | (grepl("nMigration_", parameter) & model == "MASCOT-GLM")) %>%
    mutate(parameter = gsub('.*nMigration_', '', parameter)) %>%
    left_join(., spatialCounts, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol"="protocol", 
                                       "parameter" ="parameterChar")) %>%
    mutate(value.y = ifelse(is.na(value.y), 0, value.y)) %>%
    select(all_of(colNames), X50, value.y, parameter, nSim, nSeq, protocol, model) %>%
    dplyr::rename(point = value.y, smooth_value = X50, upper = colNames[1], lower = colNames[2]) %>%
    mutate(ir = x)
})

# Plots
n_r = ifelse(grepl("7", runType), 7, 3)
spatialCounts %>%
  filter(!parameterChar %in% paste0("Region", 1:n_r, "_Region", 1:n_r)) %>%
  group_by(matrix, nSeq) %>%
  summarise(n_mis = sum(is.na(value) | value == 0), n_tot = n())

allPlots = list()
i = 1 
s = c("Kendall's tau", "Calibration", "Mean relative\n95% HPD width", 
      "Mean relative bias", "Weighted Interval score")

for (statistic in s) {
  
  if (statistic == "Kendall's tau") statistic_df = rename(regressionMJ, stat = tau)
  if (statistic == "Calibration") statistic_df = rename(migrationEventsCalibration, stat = perc)
  if (statistic == "Mean relative\n95% HPD width") statistic_df = migrationEvents %>%
      filter(value.y != 0) %>%
      group_by(protocol, nSeq, model) %>%
      summarise(stat = mean( (X97_5_hpd - X2_5_hpd) / value.y, na.rm = T))
  if (statistic == "Mean relative bias") statistic_df = migrationEvents %>%
      filter(value.y != 0) %>%
      group_by(protocol, nSeq, model) %>%
      summarise(stat = mean( (X50 - value.y) / value.y, na.rm = T))
  
  if (statistic == "Weighted Interval score") {
    statistic_df = data %>%
      filter((grepl("bssvs_nMigration_", parameter) & BF >=3) | (grepl("nMigration_", parameter) & model == "MASCOT-GLM")) %>%
      mutate(parameter = gsub('.*nMigration_', '', parameter)) %>%
      left_join(., spatialCounts, by = c("nSim" = "nSim", "nSeq" = "nSeq", "protocol"="protocol", 
                                         "parameter" ="parameterChar")) %>%
      mutate(value.y = ifelse(is.na(value.y), 0, value.y)) %>%
      mutate(stat = WIS(wis_data), nSeq = recode(nSeq, `150` = "150 sequences", `500` = "500 sequences")) 
    
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
                                 labels = c("uni. surv.", "region", "region+year"))) %>%
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
        geom_line(size = 0.8, alpha = 0.7) +
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
          geom_line(size = 0.8, alpha = 0.7) +
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
  
  if (statistic == "Kendall's tau") {
    ymin = ifelse(min(statistic_df$stat) < 0, min(statistic_df$stat), 0)
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(ymin, 1)
  }
  if (statistic == "Calibration") {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(c(0,100))
  }
  
  if (statistic %in% c("Mean relative bias", "Mean relative\n95% HPD width")) {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + ylim(min(statistic_df$stat),max(statistic_df$stat))
  }
  if (!statistic %in% c("Weighted Interval score", "Mean relative bias")) {
    for (l in i:(i+2)) allPlots[[l]] = allPlots[[l]] + rremove("x.text") + rremove("x.ticks")
  }
  
  i = i+3
}

# All plots
p1 = ggarrange(plotlist = allPlots[c(1:12)], ncol = 3, nrow = 4, 
               labels = c("A", "F", "", "B", "G", "", "C", "H", "", "D", "I", ""), 
               widths = rep(c(1.3,1,1), 4), legend = "top", common.legend = T,
               align = "hv")

p2 = ggarrange(plotlist = allPlots[c(13:15)], ncol = 3,
               labels = c("E", "J", ""), widths = c(1.5,1,1), 
               vjust = 1, hjust = 0, align = "hv", legend = "none")

tosave = ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(3.3,1), align = "hv", common.legend = T, legend = "bottom")

ggsave(paste0("2.Figures/", runType, "/total_migration_counts.png"), tosave, width = 9, height = 12)
