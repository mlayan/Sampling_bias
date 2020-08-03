#############################################################
##          ANALYSIS OF BEAST1 AND BEAST2 RUNS
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

#############################################################
## Load data 
## Selected runs 
selected = read.table("2.Figures/2.selected_runs.txt", stringsAsFactors = FALSE, 
                      header = TRUE, sep = "\t")

## Tree topology regression data
topologies1 = read.table('HKY_M1/analyses/topology_regression_sim1to10.txt', 
                         stringsAsFactors = FALSE, header = TRUE, sep = "\t")
topologies2 = read.table('HKY_M2/analyses/topology_regression_sim1to10.txt', 
                         stringsAsFactors = FALSE, header = TRUE, sep = "\t")
topologies = bind_rows(topologies1, topologies2) %>%
  mutate(model = recode(model, `mascot_v7` = "mascot")) %>%
  right_join(., selected, by = c("model", "nSim", "protocol", "nSeq", "matrix")) %>%
  mutate(protocol = factor(protocol, c(bias, surveillance)))

         
#############################################################
## Plots
# Surveillance
ggplot(topologies[topologies$protocol %in% surveillance, ], aes(x = protocol, y = coef, fill = model)) +
  facet_grid(nSeq ~ matrix, labeller = "label_both") +
  geom_boxplot() +
    scale_fill_manual(values = c(wes_palette("Darjeeling1")[2:3])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sampling protocol", y = "Determination coefficient R2", col = "BEAST model") 
ggsave("2.Figures/tree_topology_surveillance.png", width = 5, height = 5)

# Bias
ggplot(topologies[topologies$protocol %in% bias, ], aes(x = protocol, y = coef, fill = model)) +
  facet_grid(nSeq ~ matrix, labeller = "label_both") +
  geom_boxplot() +
  scale_fill_manual(values = c(wes_palette("Darjeeling1")[2:3])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sampling protocol", y = "Determinaton coefficient R2", col = "BEAST model") 
ggsave("2.Figures/tree_topology_bias.png", width = 7, height = 5)

############################################################# 
## Paired wilcoxon
topologies_wide = select(topologies, model, coef, nSim, protocol, nSeq, matrix) %>%
  spread(., model, coef)

args = expand.grid(c(150, 500), 1:2, c(bias, surveillance))
p_values = data.frame(t(apply(args, 1, tests)))
colnames(p_values) = c("nSeq", "matrix", "protocol", "p_value")
p_values$p_value = as.numeric(p_values$p_value)
p_values[p_values$p_value < 0.05,]
  
 