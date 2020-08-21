#############################################################
##            FUNCTION FOR SAMPLING PROTOCOLS
##                        
##
## Creation Date: 11-02-2020
## Last Update: 15-05-2020
## Maylis Layan
#############################################################

# Surveillance samples
sampling_protocols <- function(metadata, surveillanceSize, surveillanceBias, 
                               sampleSizes, protocols, lowerYear, index.case, 
                               directory, nSim) {
  
  ######################################################
  ## Prepare metadata
  # Add a column for the surveillance bias
  metadata$surveillanceBias <- sapply(metadata$patch, function(x) surveillanceBias[x, 2])
  
  # Add a year column
  metadata$year <- floor(metadata$samp.dates)
  
  ######################################################
  # Surveillance sample
  initRegion <- metadata$patch[metadata$short.names %in% index.case]
  surveillance <- sample(metadata$short.names[metadata$year >= lowerYear],
                         surveillanceSize, 
                         prob = metadata$surveillanceBias[metadata$year >= lowerYear])
  metadata$surveillance <- as.numeric(metadata$short.names %in% surveillance)
    
    # Regions and years represented by the surveillance sample
    rSurveillance <- sort(unique(metadata$patch[metadata$surveillance == 1]))
    nRegions <- length(rSurveillance)
    ySurveillance <- sort(unique(metadata$year[metadata$surveillance == 1]))
    nYears <- length(ySurveillance)
    
  # List of sequences sampled by the surveillance system per region
  seqR <- lapply(rSurveillance, function(x) metadata$short.names[metadata$surveillance == 1 & metadata$patch == x])
  
  # List of sequences sampled by the surveillance system per region per year
  seqRY <- vector("list", nRegions*nYears)
  names(seqRY) <- apply(expand.grid(rSurveillance, ySurveillance), 1, paste, collapse = "_")
  
  for (r in rSurveillance) {
    for (y in ySurveillance) {
      seqRY[[paste0(r,"_", y)]] <- metadata$short.names[metadata$patch == r & metadata$year == y & metadata$surveillance == 1]
    }
  }
  
  ######################################################
  # Samples
  if ("all" %in% protocols) protocols <- c("uniform", "stratified", "maxPerRegion", "maxPerRegionYear")
  
  # Uniform sample - force location representativness
  if ("uniform" %in% protocols) {
    
    for (s in sampleSizes) {
      col.name <- paste0("uniformS_", s)
      metadata[[col.name]] <- NA
      uniform <- c()
      while (length(table(metadata$patch[metadata[[col.name]] == 1])) < nRegions) {
        uniform <- sample(surveillance, s)
        metadata[[col.name]] <- as.numeric(metadata$short.names %in% uniform)
      } 
    }

  }

  
  # Stratified sample
  if ("stratified" %in% protocols) {
    
    freqTab <- table(metadata$patch[metadata$surveillance == 1]) / surveillanceSize
    metadata$freq <- sapply(metadata$patch, function(x) freqTab[as.character(x)])
    metadata$freq[metadata$surveillance == 0] <- 0
    
    for (s in sampleSizes) {
      col.name <- paste0("stratified_", s)
      stratified <- c()
      stratified <- sample(metadata$short.names, s, prob = metadata$freq)
      metadata[[col.name]] <- as.numeric(metadata$short.names %in% stratified)
    }

  }
  
  
  # Maximum per region
  if ("maxPerRegion" %in% protocols) {
    
    for (s in sampleSizes) {
      col.name <- paste0("maxPerRegion_", s)
      
      # Number of sequences to draw from each patch
      perRegion <- floor(s / nRegions)
      left <- s %% nRegions
      
      # Get all samples from under-represented locations
      # and sample perRegion sequences from over-represented regions
      maxPerRegion <- NULL
      maxPerRegion <- lapply(seqR, function(x, N) 
        if(length(x) > N) {
          sample(x, size = perRegion, replace = FALSE)
        } else {
          x
        }, N = perRegion)

      # Number of sequences left to sample
      supp <- s - length(unlist(maxPerRegion))

      # Draw supp for regions that are more represented
      available <- mapply(function(x,y) x[!x %in% y], seqR, maxPerRegion)
      while(supp > 0) {

        # Regions with sequences left
        regionsThatCanBeSampled <- as.numeric(sapply(available, length) > 0)
        regionToBeSampled <- sample(rSurveillance, size = 1, prob = regionsThatCanBeSampled)

        # Sample
        newSample <- ifelse(length(available[[regionToBeSampled]]) > 1,
                            sample(available[[regionToBeSampled]], 1),
                            available[[regionToBeSampled]])

        # Push back the list
        maxPerRegion[[regionToBeSampled]] <- c(maxPerRegion[[regionToBeSampled]], newSample)
        
        # Update available
        available[[regionToBeSampled]] <- available[[regionToBeSampled]][!available[[regionToBeSampled]] %in% newSample]
        
        # Update supp
        supp <- supp - 1
      }

      # Add a maximumPerRegion column
      metadata[[col.name]] <- as.numeric(metadata$short.names %in% unlist(maxPerRegion))

    }
    
  }

  
  # Maximum per region and per year
  if ("maxPerRegionYear" %in% protocols) {
    
    for (s in sampleSizes) {
      # Column
      col.name <- paste0("maxPerRegionYear_", s)
      
      # Number of sequences to draw from each patch and each year
      perRY <- floor(s / (nRegions*nYears))
      leftRY <- s %% (nRegions*nYears)
      
      # Get all samples from under-represented locations
      # and sample perRegion sequences from over-represented regions
      maxPerRegionYear <- NULL
      maxPerRegionYear <- sapply(seqRY, function(x) if(length(x) > perRY) {sample(x, perRY, replace = FALSE)} else {x})
      
      # Number of sequences left to sample
      suppRY <- s - length(unlist(maxPerRegionYear))
      
      # Draw supp for regions that are more represented
      availableRY <- mapply(function(x,y) x[!x %in% y], seqRY, maxPerRegionYear)
      while(suppRY > 0) {
        
        # Regions with sequences left
        regionYearAvailable <- (sapply(availableRY, length) > 0)
        regionYearToBeSampled <- sample(names(regionYearAvailable), size = 1, prob = as.numeric(regionYearAvailable))
        
        # Sample
        newSample <- ifelse(length(availableRY[[regionYearToBeSampled]]) > 1, 
                            sample(availableRY[[regionYearToBeSampled]], 1), 
                            availableRY[[regionYearToBeSampled]])
        
        # Push back the list
        maxPerRegionYear[[regionYearToBeSampled]] <- c(maxPerRegionYear[[regionYearToBeSampled]], newSample)
        
        # Update availableRY
        availableRY[[regionYearToBeSampled]] <- availableRY[[regionYearToBeSampled]][!availableRY[[regionYearToBeSampled]] %in% newSample]
        
        # Update suppRY
        suppRY <- suppRY - 1
        
      }
      
      # Add a maximumPerRegion column
      metadata[[col.name]] <- as.numeric(metadata$short.names %in% unlist(maxPerRegionYear))
    }
   
  }
  
  # Write summary table
  coverage <- as.matrix(expand.grid(patch = rSurveillance, year = ySurveillance))
  colNames <- apply(expand.grid(protocols, sampleSizes), 1, paste, collapse = "_")
  summaryDF <- apply(coverage, 1, function(x)
    sapply(colNames, 
           function(y) nrow(metadata[metadata[[y]] == 1 & metadata$patch == x[1] & metadata$year == x[2], ])))
  
  summaryDF <- data.frame(cbind(t(summaryDF), coverage))
  write.table(summaryDF, paste0(directory, "/files/sim", nSim, "_sample_composition.txt"), row.names = FALSE, sep ="\t")
  
  # Return updated metadata data frame
  return(metadata)
}




# Systematic biased sampling
sample.seq <- function(samp.size, w, index.case, simulationNb, 
                       snp, metadata, directory) {
  
  # Draw samples
  if (max(w) == 1) {
    
    samp.type <- "uniform"
    file.name <- paste(samp.type, samp.size, sep = "_")
    selection.index <- c()
    
    while (!1 %in% metadata$patch[selection.index]) {
      selection.index <- sample(1:length(w), samp.size) 
    }
    
  } else {
    
    samp.type <- "biased"
    w.name <- max(w)
    file.name <- paste(samp.type, w.name, samp.size, sep = "_")
    selection.index <- c()
    
    while (!1 %in% metadata$patch[selection.index]) {
      selection.index <- sample(1:length(w), samp.size, prob = w)   
    }
    
  }
  
  # Fasta file
  selection.snp <- snp[selection.index]
  selection.sequences <- lapply(selection.snp, write.sample.seq, index.case = index.case, snp = snp)
  selection.names <- metadata$traits[selection.index]
  write.fasta(sequences = selection.sequences, names = selection.names, 
              file.out = paste0(directory, "/files/sim", simulationNb, "_sequences_", file.name, ".fasta"))
  
  # Geographic traits file
  traits <- metadata[selection.index, c("traits", "regions")]
  write.table(traits, paste0(directory, "/files/sim", simulationNb, "_traits_", file.name, ".txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Return vector of selected sequences
  out <- rep(0, nrow(metadata))
  out[selection.index] <- 1
  return(out)
}


##########################################
## Write sequences in sample.seq function
write.sample.seq <- function(x, index.case, snp) {
  out <- snp[[index.case]]$nuc
  out[x$loc] <- x$nuc
  return(out)
}


