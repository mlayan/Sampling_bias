#############################################################
##        SIMULATION OF REALISTIC EPIDEMICS IN
##                        MOROCCO
##
##                    MAIN FUNCTION
##
## Creation Date: 06-01-2020
## Last Update: 10-11-2020
## Maylis Layan
#############################################################

rm(list = ls())
setwd("7demes")

# Arguments 
# Give a value between 1 and 5 to determine which seed should be used 
args <- commandArgs(trailingOnly=T)
if (length(args) == 0) {
  stop("Need arguments!")
}

# Print the % of migration
print(paste0("Simulations with ", args[2], "% of migration"))
migRates <- c(0.9, 1.7, 2.6, 3.6, 5.5)
names(migRates) <- c(0.5, 1, 1.5, 2, 3)

# Librairies
library(dplyr)
library(Rcpp)
library(seqinr)
library(doParallel)
library(doRNG)
library(foreach)

# Models and functions
# Epidemic model
sourceCpp("../R/cpp_model_v5.cpp")
sourceCpp("../R/cpp_beta_mat.cpp")

# Evolutionary model
source('../R/HKY_v2.R')

# Sampling protocols
source('../R/sampling_protocols_v2.R')

# Arguments
  # Start of the epidemic 
initPatch <- 1

  # Sequence simulation
ref <- read.fasta("../inputfiles/RV2627.fasta", 
	forceDNAtolower = FALSE, 
	set.attributes = FALSE)[[1]]
mut.rate <- 2.44e-4/365.25    # See Troupin et al., 2016
kappa <- 2
root.date <- as.Date("1989-01-01")

  # Sampling schemes
regions = c("Region1", 
            "Region2", 
            "Region3",
            "Region4",
            "Region5",
            "Region6",
            "Region7")
sampleSizes <- c(150, 500)
lowerYear = 1990
biasWeights = c(2.5, 5, 10, 20, 50)

	# Surveillance bias
surveillanceSize <- 5000
surveillanceProtocols = c("uniform", "maxPerRegion", "maxPerRegionYear") 

	# Systematic bias
weightsRegions <- matrix(1, ncol = length(biasWeights) + 1, nrow = length(regions))
weightsRegions[3, 2:6] <- biasWeights # Region3 will be oversampled
weightsRegions[4, 2:6] <- biasWeights # Region4 will be oversampled

rownames(weightsRegions) <- regions

# Human population
H <- read.csv("../inputfiles/inhabitants_per_region_UNadj.csv", 
	stringsAsFactors = FALSE, 
	fileEncoding="latin1")
H$aggregation[H$region %in% c("Tanger - Tétouan",
                              "Taza - Al Hoceima - Taounate")] <- "Region1"

H$aggregation[H$region == "Oriental"] <- "Region2"

H$aggregation[H$region %in% c("Gharb - Chrarda - Béni Hssen",
                              "Rabat - Salé - Zemmour - Zaer",
                              "Grand Casablanca",
                              "Chaouia - Ouardigha")] <- "Region3"

H$aggregation[H$region %in% c("Doukkala - Abda",
                                 "Marrakech - Tensift - Al Haouz")] <- "Region4"

H$aggregation[H$region %in% c("Meknès - Tafilalet",
                              "Fès - Boulemane")] <- "Region5"

H$aggregation[H$region %in% c("Souss - Massa - Draâ",
                              "Tadla - Azilal")] <- "Region6"

H$aggregation[H$region %in% c("Guelmim - Es-Semara",
                              "Laâyoune - Boujdour - Sakia El Hamra")] <- "Region7"

human_pop = group_by(H, aggregation) %>% summarise(pop = sum(ind)) %>% .$pop

  # Human mobility
distance = as.matrix(read.table("../inputfiles/mobility_matrix_ecoregions.txt", 
	header = FALSE))
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
connection_strength <- migRates[paste0(args[2])] * 1e-8
time.max <- round(365.25 * n.years, 0)
proba <- dgamma(1:15, shape = 3, scale = 1.1)
inf_period_dist <- proba/sum(proba)

# beta_mat = cpp_beta_adj(human_pop, beta, distance,
#              nu.source, mu.destination, epsilon,
#              iota, connection_strength, mobility)
# 
# beta_mat <- t(t(beta_mat)*human_pop*dog.ratio)*3.3
# Ri <- apply(beta_mat, 1, sum)
# Rii <- diag(beta_mat)/Ri*100
# X <- data.frame("Rii.Ri" = Rii, Rii = diag(beta_mat), Ri = Ri)
# tail(rbind(X, c(mean(Rii), mean(diag(beta_mat)), mean(Ri))))

#-----------------------------------------------------------
# Simulation
#-----------------------------------------------------------
registerDoParallel(cores = 10)

# Set the range of the simulation index
seed <- as.integer(Sys.time() + Sys.getpid())
if (args[1] == 1) seed = 1605115739
if (args[1] == 2) seed = 1605122202
if (args[1] == 3) seed = 1605107965
if (args[1] == 4) seed = 1605118907
if (args[1] == 5) seed = 1605118168
print(paste('Random seed =', seed))
set.seed(seed)

# Set the number of simulation
print(paste("Iteration number =", as.numeric(args[1])))

lowerBound = as.numeric(args[1])*10-9
upperBound = as.numeric(args[1])*10


# Simulation
out <- foreach (i=lowerBound:upperBound,
                .combine = rbind,
                .packages = c("seqinr", "lubridate", "Rcpp", "dplyr", "ape"))  %dorng% {
                  
                  # Create directory if not existing
                  listDir <- list.dirs(full.names = FALSE)
                  directory <- paste0("simulation", i, "/files")
                  if (!directory %in% listDir) dir.create(directory, recursive = T)
                  
                  ## Epidemic simulation--------------------------------------
                  start <- Sys.time()
                  res <- cpp_model(human_pop, dog.ratio, beta, distance, init.inf,
                                   time.max, inf_period_dist, birth.rate, death.rate,
                                   nu.source, mu.destination, epsilon, iota, connection_strength, 
                                   mobility)
                  
                  while (res$size < min.cases | res$duration != res$duration) {
                    res <- cpp_model(human_pop, dog.ratio, beta, distance, init.inf,
                                     time.max, inf_period_dist, birth.rate, death.rate,
                                     nu.source, mu.destination, epsilon, iota, connection_strength, 
                                     mobility)
                    
                  }
                  end <- Sys.time()
                  print(paste("Found - simulation", i)) 
                  print(paste("Time spent", difftime(end, start, "secs"), "seconds - simulation", i))
                  
                  ## Sequence simulation---------------------------------------
                  # Simulate sequences
                  start2 <- Sys.time()
                  snp <- HKY_model(res, ref, mut.rate, kappa, root.date, directory, simulationNb = i)
                  end2 <- Sys.time()
                  print(paste("SNP done - simulation", i, " - Time spent:", difftime(end2, start2, "hours"), "hours")) 
                  
                  # Metadata table
                  metadata <- read.delim(paste0(directory, "/sim", i, "_metadata.txt"), stringsAsFactors = FALSE)
                  
                  # Add names of regions
                  metadata$regions <- sapply(metadata$patch, function(x) regions[x])

                  ## Sample sequences-----------------------------------------
                  # Index case
                  index.case <- paste0(initPatch-1, "_0")
                  
                  # Sample according to different sampling protocols
                  metadata <- sampling_protocols(metadata, lowerYear, initPatch, 
                  	directory, sampleSizes, i, weightsRegions, surveillanceSize, 
                  	weightsRegions[, -1], surveillanceProtocols)
                  
                  ################################################################
                  ## DUE TO AN OVERSIGHT, IT WAS POSSIBLE TO SAMPLE SEQUENCES
                  ## DURING THE 1st YEAR OF THE EPIDEMIC IN THE SYSTEMATIC BIAS 
                  ## SCHEMES 
                  ## TO AVOID RE-RUNNING THE FULL ANALYSIS AND TO ENSURE REPRODUCIBLE 
                  ## RESULTS, WE DECIDED TO KEEP THE SAME SEED AND WE ADD THE FOLLOWING 
                  ## LINES OF CODE TO RE-SAMPLE CORRECTLY THE SEQUENCES  
                  
                  # Identify conditions
                  conditions <- colnames(metadata)[grep("_[0-9]+", colnames(metadata))]
                  conditionsToResample = c() 
                  for (cond in conditions) {
                    if(any(metadata$samp.dates[metadata[[cond]] == 1] < lowerYear)) conditionsToResample = c(conditionsToResample, cond) 
                  }
                  
                  # Index case
                  indexCase <- paste0(initPatch-1, "_0")
                  
                  # Systematic bias
                  # Names of biased protocols
                  biased <- paste0("biased_", apply(weightsRegions, 2, max))
                  biased[biased == "biased_1"] <- "uniform"
                  colnames(weightsRegions) <- biased
                  
                  # Draw samples
                  for (cond in conditionsToResample) {
                    COL = gsub("_[0-9]+$", "", cond)
                    s = as.numeric(gsub("^.*_", "", cond))
                    weightsCOL <- sapply(metadata$patch, function(x) weightsRegions[x, COL])
                    weightsCOL[is.na(weightsCOL)] = 0
                    selectedIndices <- c()
                    
                    # Draw samples by forcing the origin of the epidemic to be sampled 
                    while (!initPatch %in% metadata$patch[selectedIndices] | any(metadata$samp.dates[selectedIndices] < lowerYear)) {
                      selectedIndices <- sample(1:nrow(metadata), s, prob = weightsCOL)
                    }
                    
                    # Add columns to metadata
                    metadata[[cond]] <- 0
                    metadata[selectedIndices, cond] <- 1
                  }
                  
                  ## THE INITIAL CODE RESUMES HERE
                  ################################################################
                  
                  # Write sequences
                  conditions <- colnames(metadata)[grep("_[0-9]+", colnames(metadata))]
                  for (cond in conditions){
                    # Sequences
                    sequences <- lapply(metadata$short.names[metadata[[cond]] == 1], 
                    	write.seq, 
                    	index.case = index.case, 
                    	snp = snp)
                    
                    write.fasta(sequences, 
                    	metadata$traits[metadata[[cond]] == 1],
                    	paste0(directory, "/sim", i, "_sequences_", cond, ".fasta"))

                    # Geographic traits
                    traits <- metadata[metadata[[cond]] == 1, c("traits", "regions")]
                    write.table(traits, 
                    	paste0(directory, "/sim", i, "_traits_", cond, ".txt"),
                    	sep = "\t", 
                    	quote = FALSE, 
                    	row.names = FALSE)
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
                  write.table(transmission.chain, paste0(directory, "/sim", i, "_transmission_chain.txt"), quote = FALSE, sep = "\t")
                  
                  
                  ## Write the description of spatial sampling--------------------------------
                  transmission.chain %>%
                    dplyr::select(patch, contains("_")) %>%
                    group_by(patch) %>%
                    summarise_all(sum, na.rm = TRUE) %>%
                    as.data.frame() ->
                    spatial.sampling

                  addedRows <- nrow(spatial.sampling) - n.patch
                  spatial.sampling$regions <- c(regions, rep(NA, addedRows))
                  
                  write.table(spatial.sampling, paste0(directory, "/sim", i, "_spatial_sampling.txt"), quote = FALSE, sep = "\t")

                  # Output: size of the epidemic
                  c(i, res$size, as.numeric(difftime(end, start, "secs")), as.numeric(difftime(end2, start2, "hours")))
                  

}

simulationsSummary <- data.frame(out)
colnames(simulationsSummary) <- c("simulation", "epidemicSize", "epidemicCompTime", "evolutionCompTime")

if (!dir.exists("generate")) dir.create("generate", showWarnings = F)
write.table(simulationsSummary, paste0("generate/simulation_summary_launch_", args[1], ".txt"), sep ="\t", row.names = FALSE)

