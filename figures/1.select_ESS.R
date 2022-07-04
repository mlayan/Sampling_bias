#############################################################
##          SELECT RUNS USING A THRESHOLD ON ESS
##
## Maylis Layan
#############################################################

rm(list = ls())
library(tidyverse)

setwd("X:/")
source("R_Functions/plot_results.R")

# ESS threshold on prior, posterior and likelihood
threshold = 200

# Simulation and inference frameworks
runTypes = c("dta", "basta", "mascot")
category = c("3demes", "7demes")
subcategory = as.character(seq(1,50,10))
nSeq = as.character(c(150,500))

# Create output directory if it doesn't exists
for (directory in category) {
  if (!dir.exists(paste0("2.Figures/", directory))) 
    dir.create(paste0("2.Figures/", directory), recursive = T)
}

#############################################################
## Load data 
#############################################################
## BEAST data
data = data.frame()

# Directories
dirs = paste0("HKY_", category, "/analyses/")

# File names
allConds = expand.grid(runTypes, nSeq, subcategory)
allConds = apply(allConds, 1, paste, collapse = "_")
allCondsGLM = expand.grid(nSeq, subcategory)
allCondsGLM = apply(allCondsGLM, 1, paste, collapse = "_")
files = c(paste0("log_files_", allConds, ".txt"), paste0("glm_files_", allCondsGLM, ".txt"))

# Load data
for (i in dirs) {
  for (f in files) {
    fileName = paste0(i, f)
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
# Select data
#############################################################
for (directory in category) {
  outputDirectory = paste0("2.Figures/", directory, "/")
  
  # Analyze ESS values
  low_ESS = data %>%
    filter(!(parameter %in% "ageRoot" & model %in% "basta")) %>%
    filter(matrix %in% directory, !grepl("indicators|nMigration|sumNonZeroRates|rootLocation|allTransitions|KL|sumMigPredNonZero|sumNePredNonZero", parameter)) %>%
    filter((grepl("bssvs_", parameter) & BF >= 3 & round(ESS) < threshold) | (!grepl("bssvs_", parameter) & round(ESS) < threshold)) %>%
    select(matrix, nSim, nSeq, protocol, model)
  low_ESS = low_ESS[!duplicated(low_ESS), ]
  rownames(low_ESS) = NULL
  
  print(nrow(low_ESS))
  print(table(low_ESS$model, low_ESS$nSeq))
  
  low_ESS %>%
    arrange(model) %>%
    write.table(.,
                paste0(outputDirectory, "low_ESS_runs.txt"),
                sep = "\t",
                row.names = F)

  # Remove runs with ESS < threshold for the posterior, likelihood and prior parameters
  high_ESS = data %>%
    select(matrix, nSim, nSeq, protocol, model) %>%
    filter(matrix %in% directory)
  high_ESS = high_ESS[!duplicated(high_ESS), ]
  high_ESS = anti_join(high_ESS, low_ESS)

  high_ESS %>%
    filter(!(model %in% "basta" & nSeq == 500)) %>%
    write.table(., paste0(outputDirectory, '2.selected_runs.txt'), sep = "\t",
              row.names = FALSE, col.names = TRUE)

  data %>%
    filter(matrix %in% directory) %>%
    right_join(., high_ESS) %>%
    filter(ESS >= threshold | is.na(ESS)) %>%
    data.frame() %>%
    write.table(., paste0(outputDirectory, "/2.selected_data.txt"), sep = "\t",
                row.names = FALSE, col.names = TRUE)
  
}

