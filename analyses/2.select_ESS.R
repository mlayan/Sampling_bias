#############################################################
##          SELECT RUNS USING A THRESHOLD ON ESS
##
## Creation Date: 18-06-2020
## Last Update: -06-2020
## Maylis Layan
#############################################################

rm(list = ls())
library(tidyverse)

source("R/plot_results.R")

runTypes = c("dta", "mascot_v7")
runType = "mascot_V7"
threshold = 100

#############################################################
## Load data 
## BEAST data
data = data.frame()
dirs = #"HKY_radiation" 
  c('HKY_M1', 'HKY_M2')
files = paste0("log_files_", runTypes, "_ESS.txt")

for (i in dirs) {
  for (f in files) {
    print(f)
    if (nrow(data) == 0) {
      data = read.table(paste0(i, '/analyses/', f), stringsAsFactors = FALSE, 
                        header = TRUE, sep = "\t")
    } else {
      tempdata = read.table(paste0(i, '/analyses/', f), stringsAsFactors = FALSE, 
                            header = TRUE, sep = "\t")
      data = bind_rows(data, tempdata)
      rm(tempdata)
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
  write.table(., paste0("2.Figures/", runType, "/2.selected_runs_", runType, ".txt"), sep = "\t", 
              row.names = FALSE, col.names = TRUE)

# Modify mascot_v* into mascot and  
# remove parameters with ESS < 100
right_join(data, high_ESS) %>%
  mutate(model = recode(model, runTypes='mascot')) %>%
  filter(ESS >= threshold || is.na(ESS)) %>%
  data.frame() %>%
  write.table(., paste0("2.Figures/", runType, "/2.selected_data_", runType, ".txt"), sep = "\t", 
              row.names = FALSE, col.names = TRUE)
  
