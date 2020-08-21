#############################################################
##          SELECT RUNS USING A THRESHOLD ON ESS
##
## Creation Date: 18-06-2020
## Last Update: -06-2020
## Maylis Layan
#############################################################

rm(list = ls())
library(tidyverse)

setwd("/mnt/gaia/MMMI_Rage/")
source("R_Functions/plot_results.R")

runTypes = c("dta", "mascot_v8")
category = "radiation_bis"
if (!dir.exists(paste0("2.Figures/", category))) dir.create(paste0("2.Figures/", category))
threshold = 200

#############################################################
## Load data 
## BEAST data
data = data.frame()
dirs = paste0("HKY_", category)
files = paste0("log_files_", runTypes, "_ESS.txt")

for (i in dirs) {
  for (f in files) {
    fileName = paste0(i, '/analyses/', f)
    if (file.exists(fileName)) {
      print(f)
      dataTemp = read.table(fileName, stringsAsFactors = FALSE,
                            header = TRUE, sep = "\t")
      data = bind_rows(data, dataTemp)
      rm(dataTemp)
    }
  }
}

#############################################################
# Analyze ESS values
low_ESS = data[data$ESS < threshold & data$parameter %in% c("likelihood", "posterior", "prior"), 
               c("matrix", "nSim", "nSeq", "protocol", "model")]
low_ESS = low_ESS[!duplicated(low_ESS), ]
rownames(low_ESS) = NULL
nrow(low_ESS)
low_ESS

# Remove runs with ESS < threshold for the posterior, likelihood and prior parameters
high_ESS = data[,c("matrix", "nSim", "nSeq", "protocol", "model")]
high_ESS = high_ESS[!duplicated(high_ESS), ]
high_ESS = anti_join(high_ESS, low_ESS)

mutate(high_ESS, model = recode(model, runTypes='mascot')) %>%
  write.table(., paste0("2.Figures/", category, "/2.selected_runs.txt"), sep = "\t", 
              row.names = FALSE, col.names = TRUE)

# Modify mascot_v* into mascot and  
# remove parameters with ESS < 100
right_join(data, high_ESS) %>%
  mutate(model = recode(model, runTypes='mascot')) %>%
  filter(ESS >= threshold || is.na(ESS)) %>%
  data.frame() %>%
  write.table(., paste0("2.Figures/", category, "/2.selected_data.txt"), sep = "\t", 
              row.names = FALSE, col.names = TRUE)
  
