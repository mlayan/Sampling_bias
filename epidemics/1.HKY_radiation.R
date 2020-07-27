#############################################################
##        SIMULATION OF REALISTIC EPIDEMICS IN
##                        MOROCCO
##
##                    MAIN FUNCTION
##
## Creation Date: 06-01-2020
## Last Update: -01-2020
## Maylis Layan
#############################################################

rm(list = ls())
#setwd("/mnt/gaia/MMMI_Rage/HKY_radiation_bis/")
setwd("..")

args <- commandArgs(trailingOnly=T)
if (length(args) == 0) {
  stop("Need arguments!")
}

# Librairies
library(dplyr)
library(Rcpp)
library(seqinr)
library(doParallel)
library(doRNG)
library(foreach)

# Models and functions
# Epidemic model
sourceCpp("../R_Functions/cpp_model_v5.cpp")
sourceCpp("../R_Functions/cpp_beta_mat.cpp")

# Evolutionary model
source('../R_Functions/HKY_v2.R')

# Sampling protocols
source('../R_Functions/sampling_protocols.R')

# Arguments
  # Spatial resolution
initPatch <- 1

  # Sequence simulation
ref <- read.fasta("../Input_files/RV2627.fasta", forceDNAtolower = FALSE, set.attributes = FALSE)[[1]]
mut.rate <- 2.44e-4/365.25
kappa <- 2
root.date <- as.Date("1989-01-01")

  # Surveillance bias
regions = c("Western Mediterranean coast",
            "Eastern Mediterranean coast",
            "North Atlantic coast",
            "South Atlantic coast",
            "Northern Atlas",
            "Southern Atlas",
            "Anti Atlas")

surveillanceSize <- 5000
surveillanceBias <- matrix(1:length(regions), ncol = 2, nrow = length(regions))
surveillanceBias[,2] <- c(1,1,5,5,1,1,1)
sampleSizes <- c(150, 500)
protocols <- "all" 
lowerYear <- 1990

  # Systematic bias
weightsRegions <- matrix(1, ncol = 6, nrow = length(regions))
weightsRegions[c(3,4),2] <- 2.5
weightsRegions[c(3,4),3] <- 5
weightsRegions[c(3,4),4] <- 10
weightsRegions[c(3,4),5] <- 20
weightsRegions[c(3,4),6] <- 50
rownames(weightsRegions) <- regions

  # Human population
inhab = read.csv("../Input_files/inhabitants_per_ecoregion_UNadj.csv", 
                 stringsAsFactors = FALSE, encoding = "latin1")
human_pop = sapply(sort(unique(inhab$agg.id)), function(x) sum(inhab$ind[inhab$agg.id == x]))

  # Human mobility
distance =  as.matrix(read.table("../Input_files/mobility_matrix_ecoregions.txt", header = FALSE))
mobility = TRUE

############################################################
## Generate epidemic
#############################################################
# Parameters
# Epidemic simulation
n.patch <- length(human_pop)
init.inf <- rep(0, n.patch)
init.inf[initPatch] <- 1
n.years <- 30
underreporting <- 0.2
cases_per_year <- 400
min.cases <- cases_per_year*n.years/underreporting
dog.ratio <- 0.1
birth.rate <- 1/365.25
death.rate <- 1/365.25
beta <- 3.2/human_pop
nu.source <- 0
mu.destination <- 0
epsilon <- 1
iota <- 0
connection_strength <- 0.00000001
time.max <- round(365.25 * n.years, 0)
proba <- dgamma(1:15, shape = 3, scale = 1.1)
inf_period_dist <- proba/sum(proba)

beta_mat = cpp_beta_adj(human_pop, beta, distance,
             nu.source, mu.destination, epsilon,
             iota, connection_strength, mobility)

beta_mat <- t(t(beta_mat)*human_pop*dog.ratio)*3.3
Ri <- apply(beta_mat, 1, sum)
Rii <- diag(beta_mat)/Ri*100
X <- data.frame("Rii.Ri" = Rii, Rii = diag(beta_mat), Ri = Ri)
tail(rbind(X, c(mean(Rii), mean(diag(beta_mat)), mean(Ri))))

#-----------------------------------------------------------
# Simulation
#-----------------------------------------------------------
registerDoParallel(cores = 10)

# Set the range of the simulation index
seed <- as.integer(Sys.time() + Sys.getpid())
print(paste('Random seed =', seed))
set.seed(seed)

# Set the number of the simulation
print(paste("Iteration number =", as.numeric(args[1])))

l <- as.numeric(args[1])*10-9
r <- as.numeric(args[1])*10


# Simulation
out <- foreach (i=l:r,
                .combine = rbind,
                .packages = c("seqinr", "lubridate", "Rcpp", "dplyr", "ape"))  %dorng% {
                  
                  # Create directory if not existing
                  listDir <- list.dirs(full.names = FALSE)
                  directory <- paste0("simulation", i)
                  if (!directory %in% listDir) dir.create(paste0(directory, "/files"), recursive = T)
                  
                  ## Epidemic simulation--------------------------------------
                  start <- Sys.time()
                  res <- cpp_model(human_pop, dog.ratio, beta, distance, init.inf,
                                   time.max, inf_period_dist, birth.rate, death.rate,
                                   nu.source, mu.destination, epsilon, iota, connection_strength, 
                                   mobility)
                  
                  while (res$size < min.cases) {
                    res <- cpp_model(human_pop, dog.ratio, beta, distance, init.inf,
                                     time.max, inf_period_dist, birth.rate, death.rate,
                                     nu.source, mu.destination, epsilon, iota, connection_strength, 
                                     mobility)
                    
                  }
                  end <- Sys.time()
                  print(paste("Found - simulation", i)) 
                  print(paste("Time spent", difftime(end, start, "sec"), "seconds - simulation", i))
                  
                  ## Sequence simulation---------------------------------------
                  # Simulate sequences
                  start2 <- Sys.time()
                  snp <- HKY_model(res, ref, mut.rate, kappa, root.date, simulationNb = i)
                  end2 <- Sys.time()
                  print(paste("SNP done - simulation", i, " - Time spent:", difftime(end2, start2, "hours"), "hours")) 
                  
                  # Metadata table
                  metadata <- read.delim(paste0(directory, "/files/sim", i, "_metadata.txt"), stringsAsFactors = FALSE)
                  
                  # Add names of regions
                  metadata$regions <- sapply(metadata$patch, function(x) regions[x])

                  ## Sample sequences-----------------------------------------
                  # Index case
                  index.case <- paste0(initPatch-1, "_0")
                  
                  # Sample according to different sampling protocols
                  metadata <- sampling_protocols(metadata, surveillanceSize, surveillanceBias, 
                                                 sampleSizes, protocols, lowerYear, index.case,
                                                 directory, i)
                  
                  # Write sequences
                  conditions <- colnames(metadata)[grep("_[0-9]+", colnames(metadata))]
                  for (cond in conditions){
                    
                    nSeq <- regmatches(cond, gregexpr("[510]+$", cond))[[1]]
                    prot <- regmatches(cond, gregexpr("^[a-zA-Z]+", cond))[[1]]
                    
                    # Sequences
                    sequences <- lapply(metadata$short.names[metadata[[cond]] == 1], write.seq, index.case = index.case, snp = snp)
                    write.fasta(sequences, metadata$traits[metadata[[cond]] == 1], 
                                paste0(directory, "/files/sim", i, "_sequences_", prot, "_", nSeq, ".fasta"))
                    
                    # Geographic traits
                    traits <- metadata[metadata[[cond]] == 1, c("traits", "regions")]
                    write.table(traits, paste0(directory, "/files/sim", i, "_traits_", prot, "_", nSeq, ".txt"), 
                                sep = "\t", quote = FALSE, row.names = FALSE)
                  }
                  
                  # Systematically biased samples
                  weights <- matrix(NA, ncol = ncol(weightsRegions), nrow = nrow(metadata))
                  weights[,1] <- 1
                  
                  for (s in 2:ncol(weightsRegions)) {
                    weights[,s] <- sapply(metadata$patch, function(x) weightsRegions[x, s])
                  }
                  
                  # Write fasta and traits files
                  for (s in sampleSizes) {
                    for (j in 1:ncol(weights)) {
                      
                      if (max(weights[,j]) != 1) {
                        col.name <- paste0("biased_", max(weights[,j]), "_", s)
                      } else{
                        col.name <- paste0("uniform_", s)
                      }
                      
                      metadata[[col.name]] <- sample.seq(s, weights[,j], index.case, i, snp, metadata, directory)
                      
                    }
                  }
                  

                  ## Write the transmission chain--------------------------------
                  # Convert epidemic output to data frame
                  transmission.chain <- data.frame(do.call("rbind", res$sourcing$data))
                  names(transmission.chain) <- c("patch.source", "ID.source", "infection.time", "death.time",
                                                 "infectious.time", "infectious.period", "generation.time", "dead", "hardstop")

                  # Add the case and the source for each infection event
                  transmission.chain$case <- rownames(transmission.chain)
                  transmission.chain$source <- paste(transmission.chain$patch.source, transmission.chain$ID.source, sep = "_")
                  
                  # Modify the patch.source counting
                  transmission.chain$patch.source <- sapply(transmission.chain$patch.source, function(x) 
                    ifelse(x < 0, x, x + 1))

                  # Add the sampled cases per sampling protocol
                  transmission.chain <- left_join(transmission.chain, metadata, by = c("case" = "short.names"))
                  transmission.chain[grep("_", colnames(transmission.chain))][is.na(transmission.chain[grep("_", colnames(transmission.chain))])] <- 0
                  
                  # Write the transmission chain
                  write.table(transmission.chain, paste0(directory, "/files/sim", i, "_transmission_chain.txt"), quote = FALSE, sep = "\t")
                  
                  
                  ## Write the description of spatial sampling--------------------------------
                  transmission.chain %>%
                    dplyr::select(patch, contains("_")) %>%
                    group_by(patch) %>%
                    summarise_all(sum, na.rm = TRUE) %>%
                    as.data.frame() ->
                    spatial.sampling

                  addedRows <- nrow(spatial.sampling) - n.patch
                  spatial.sampling$regions <- c(regions, rep(NA, addedRows))
                  
                  write.table(spatial.sampling, paste0(directory, "/files/sim", i, "_spatial_sampling.txt"), quote = FALSE, sep = "\t")

                  # Output: size of the epidemic
                  c(i, res$size, as.numeric(difftime(end, start, "sec")), as.numeric(difftime(end2, start2, "hours")))
                  

}

simulationsSummary <- data.frame(out)
colnames(simulationsSummary) <- c("simulation", "epidemicSize", "epidemicCompTime", "evolutionCompTime")

write.table(simulationsSummary, paste0("simulation_summary_launch", args[1], ".txt"), sep ="\t", row.names = FALSE)

