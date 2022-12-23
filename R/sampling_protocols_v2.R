#############################################################
##            FUNCTION FOR SAMPLING PROTOCOLS
##                        
##
## Creation Date: 05-08-2020
## Last Update:
## Maylis Layan
#############################################################

# Draw samples for each sampling scheme
sampling_protocols <- function(
  metadata,
  lowerYear, 
  initPatch, 
  directory, 
  sampleSizes, 
  nSim, 
  biasedWeights = NA, 
  surveillanceSize = NA, 
  surveillanceRegions = NA, 
  surveillanceProtocols = NA
  ) 
{
  
  # Index case
  indexCase <- paste0(initPatch-1, "_0")
  
  ######################################################
  # Systematic bias
  if (!any(is.na(biasedWeights))) {
    # Names of biased protocols
    biased <- paste0("biased_", apply(biasedWeights, 2, max))
    biased[biased == "biased_1"] <- "uniform"
    colnames(biasedWeights) <- biased
    
    # Draw samples
    for (s in sampleSizes) {
      for (COL in biased) {
        weightsCOL <- sapply(metadata$patch, function(x) biasedWeights[x, COL])
        selectedIndices <- c()
        
        # Draw samples by forcing the origin of the epidemic to be sampled 
        while (!initPatch %in% metadata$patch[selectedIndices]) {
          selectedIndices <- sample(1:nrow(metadata), s, prob = weightsCOL)
        }
        
        # Add columns to metadata
        metadata[[paste0(COL, "_", s)]] <- 0
        metadata[selectedIndices, paste0(COL, "_", s)] <- 1
      }
    }
  }
  
  ######################################################
  ## Surveillance-related bias
  if (!is.na(surveillanceSize)) {
    for (n in 1:ncol(surveillanceRegions)) {

      surveillanceBias = surveillanceRegions[, n] 
      metadata = surveillanceSamples(
        metadata, 
        surveillanceSize, 
        surveillanceBias, 
        surveillanceProtocols, 
        sampleSizes, 
        initPatch, 
        lowerYear
        )

    }
  }

  # # Write summary table
  # coverage <- as.matrix(expand.grid(patch = rSurveillance, year = ySurveillance))
  # colNames <- apply(expand.grid(c(surveillanceProtocols, biased), sampleSizes), 1, paste, collapse = "_")
  # summaryDF <- apply(coverage, 1, function(x)
  #   sapply(colNames, 
  #          function(y) nrow(metadata[metadata[[y]] == 1 & metadata$patch == x[1] & metadata$year == x[2], ])))
  
  # summaryDF <- data.frame(cbind(t(summaryDF), coverage))
  # write.table(summaryDF, paste0(directory, "/sim", nSim, "_sample_composition.txt"), row.names = FALSE, sep ="\t")
  
  # Return updated metadata data frame
  return(metadata)
}








# Draw samples for surveillance-related biased schemes
surveillanceSamples <- function(
  metadata, 
  size, 
  biases, 
  protocols, 
  sampleSizes, 
  initPatch, 
  lowerYear
  ) 
{

  # Surveillance bias
  bias = max(biases)

  # Add a column for the surveillance bias
  surveillanceWeights <- sapply(metadata$patch, function(x) biases[x])

  # Add a year column
  metadata$year <- floor(metadata$samp.dates)

  ######################################################
  # Surveillance sample
  surveillanceIndices <- c()
  while (!initPatch %in% metadata[surveillanceIndices, "patch"]) {
    surveillanceIndices <- sample(rownames(metadata)[metadata$year >= lowerYear],
                                  size, 
                                  prob = surveillanceWeights[metadata$year >= lowerYear])  
  }
  metadata$surveillance <- 0
  metadata[surveillanceIndices, "surveillance"] <- 1

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
  protocols[protocols == "uniform"] <- "uniformS"
   
  # Uniform sample - force location representativeness
  if ("uniformS" %in% protocols) {
    
    for (s in sampleSizes) {
      col.name <- paste0("uniformS_", bias, "_", s)
      metadata[[col.name]] <- NA
      uniform <- c()
      while (length(table(metadata$patch[metadata[[col.name]] == 1])) < nRegions) {
        uniform <- sample(metadata$short.names[metadata$surveillance == 1], s)
        metadata[[col.name]] <- as.numeric(metadata$short.names %in% uniform)
      } 
    }
    
  }


  # # Stratified sample
  # if ("stratified" %in% protocols) {
  #   
  #   freqTab <- table(metadata$patch[metadata$surveillance == 1]) / size
  #   metadata$freq <- sapply(metadata$patch, function(x) freqTab[as.character(x)])
  #   metadata$freq[metadata$surveillance == 0] <- 0
  #   
  #   for (s in sampleSizes) {
  #     col.name <- paste0("stratified_", bias, "_", s)
  #     stratified <- c()
  #     stratified <- sample(metadata$short.names, s, prob = metadata$freq)
  #     metadata[[col.name]] <- as.numeric(metadata$short.names %in% stratified)
  #   }
  #   
  # }


  # Maximum per region
  if ("maxPerRegion" %in% protocols) {
    
    for (s in sampleSizes) {
      col.name <- paste0("maxPerRegion_", bias, "_", s)
      
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
      col.name <- paste0("maxPerRegionYear_", bias, "_", s)
      
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

  return(metadata)
}
