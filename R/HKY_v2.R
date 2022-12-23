#############################################################
##      SIMULATE SEQUENCES ACCORDING TO THE HKY MODEL  
##                        
##
## Creation Date: 28-01-2020
## Last Update: 15-05-2020
## Maylis Layan
#############################################################


# FUNCTION TO WRITE SEQUENCES
write.seq <- function(case, index.case, snp, error.rate = 0) {
  
  # Write full sequence based on reference sequence
  out <- snp[[index.case]]$nuc
  out[snp[[case]]$loc] <- unlist(snp[[case]]$nuc)
  
  if (error.rate > 1 | error.rate < 0) {
    stop("Error rate should range between 0 and 1")
  
  } else if (error.rate!= 0) {
    bases <- c("A", "T", "C", "G")
    
    # Location of errors
    errorIndex <- sample(length(out), size = floor(length(out)*error.rate), replace = FALSE)
    
    # Edit sequence with sampled error 
    out[errorIndex] <- sapply(out[errorIndex], function(x) sample(bases[!bases %in% x], 1))
    
  }
  
  return(out)
}



# HKY SUBSTITUTION MODEL
HKY_model <- function(sim, ref, mu, kappa, root.date, directory, simulationNb = NULL, coords = NULL) {
  
  #################################################
  ## PACKAGES AND PARAMETERS
  # Load packages
  require(seqinr, quietly = TRUE)
  require(lubridate, quietly = TRUE)

  # Convert cases.info from a list to a data frame
  cases.info <- data.frame(do.call("rbind", sim$sourcing$data))
  names(cases.info) <- c("patch.source", "ID.source", "infection.time", "death.time", "infectious.time", 
                       "infectious.period", "generation.time", "dead", "hard.stop")
  
  # Number of initial sequences
  n.init <- length(which(cases.info$infection.time < 0))
  names.init <- rownames(cases.info[which(cases.info$infection.time < 0),])
  
  # Informations from the sequence
  # By default, it is the first sequence appearing in cases.info
  genome.length <- length(ref)
  baseComp <- as.numeric(table(ref)) 
  baseFreq <- baseComp / genome.length
  
  # Vector of bases
  bases <- names(table(ref))
  names(baseComp) <- names(baseFreq) <- bases
  
  # Eliminate individuals who die before becoming infectious
  cases.info <- cases.info[cases.info$infectious.time >= 0, ]
  
  # Lists of SNPs 
  snp <- vector("list", nrow(cases.info))
  names(snp) <- rownames(cases.info)
  snp[[names.init[1]]] <- list(loc = c(), 
                             nuc = ref, 
                             date = decimal_date(as.Date(root.date)), 
                             nb.snp = 0, 
                             d.dist = 0)
  
  ################################################
  ## INDEX CASES
  if (n.init > 1) {
    for (i in 2:n.init) {
      loc <- sample(1:genome.length, round(0.1*genome.length, 0))
      nuc <- sapply(loc, function(x) sample(bases[!(bases %in% ref[x])],1))
      
      tempSeq <- ref
      tempSeq[loc] <- nuc
      
      snp[[names.init[i]]] <- list(loc = loc, 
                                   nuc = nuc, 
                                   date = root.date, 
                                   nb.snp = 0, 
                                   d.dist = 0,
                                   baseFreq = as.numeric(table(tempSeq)) / genome.length)
    }
  }
  
  # Time max can exceed the time.max used for the simulation if the epidemic
  # doesn't stop before the end of the simulation 
  # To speed up the evolution loop, we define only the time steps that 
  # should be treated
  time.vector <- sort(unique(cases.info$infection.time))
  time.vector <- time.vector[which(time.vector >= 0)] 
  
  ################################################
  ## EVOLUTIONARY PARAMETERS OF THE HKY MODEL
  # HKY parameters 
  beta <- mu/(2*(baseFreq['A']*baseFreq['C']+baseFreq['A']*baseFreq['T']+baseFreq['C']*baseFreq['G']+baseFreq['G']*baseFreq['T'])+
                2*kappa*(baseFreq['A']*baseFreq['G']+baseFreq['C']*baseFreq['T']))
  alpha <- kappa*beta
  
  # Q matrix
  Q <- matrix(NA, ncol = 4, nrow = 4)
  Q[upper.tri(Q)] <- c(beta*baseFreq['C'], alpha*baseFreq['G'], 
                       beta*baseFreq['G'], beta*baseFreq['T'], 
                       alpha*baseFreq['T'], beta*baseFreq['T']) 
  Q[lower.tri(Q)] <- c(beta*baseFreq['A'], alpha*baseFreq['A'], 
                       beta*baseFreq['A'], beta*baseFreq['C'], 
                       alpha*baseFreq['C'], beta*baseFreq['G'])
  diag(Q) <- -apply(Q, 1, sum, na.rm = TRUE)
   
  # Compute the invertible matrix of Q 
  # and its eigenvalues 
  E <- eigen(Q)$vectors  
  E_1 <- solve(E)
  L <- diag(eigen(Q)$values)
  
  ###############################################
  ## SEQUENCE EVOLUTION
  # Function
  evol <- function(x, t) {
    
    # Retrieve individual parameters
    vec <- cases.info[rownames(cases.info) == x, ]
    
    # Get the generation time 
    generation.time <- as.numeric(vec["generation.time"])
    
    # Compute the diagonal matrix for the 
    # current generation time
    L_curr <- L
    diag(L_curr) <- exp(diag(L)*generation.time)
    
    # Compute the transition probability matrix 
    p_t <- E %*% L_curr %*% solve(E)
    colnames(p_t) <- rownames(p_t) <- bases
    if (!all(apply(p_t, 1, sum) > rep(1-1e-10, length(bases)))) stop("Transition probabilities do not sum to 1!")
    
    # Update the composition of the sequence
    ID.infector <- paste(vec[c('patch.source', 'ID.source')], collapse = "_")
    
    if (ID.infector == names.init[1]){
      
      # Sample new sequence
      newSequence <- sapply(ref, function(x) sample(bases,
                                                    size = 1,
                                                    prob = p_t[x, ]))
      
      # Store the changes
      loc <- which(newSequence != ref)
      nuc <- newSequence[loc]
      d.dist <- length(loc)
      
    } else {
      
      # Sequence of the source
      source <- ref
      if (length(snp[[ID.infector]]$loc) > 0) source[snp[[ID.infector]]$loc] <- snp[[ID.infector]]$nuc
      
      # Sample new sequence
      newSequence <- sapply(source, function(x) sample(bases,
                                                       size = 1,
                                                       prob = p_t[x, ]))
      
      # Store the changes
      # Since simulated epidemics are large, the size of sequence alignments could easily overload the R cache. 
      # Thus, we took a similar approach as the R package SEEDY (doi: 10.1371/journal.pone.0129745) by storing 
      # SNPs that differ from the reference genome in list objects instead of storing whole genome sequences.
      # Store SNPs according to the reference genome which enables to get the SNPs of the infector sequence
      loc <- which(newSequence != ref)                
      nuc <- newSequence[loc]
      d.dist <- length(which(newSequence != source))  # Compute the number of new SNPs compared to the source sequence 
      
    }
    
    # root.date should be expressed as "yyyy-mm-dd", ie "%Y-%m-%d" or "%Y/%m/%d"
    final.date <- decimal_date(as.Date(root.date) + as.numeric(vec['infectious.time']))
    final.date <- round(final.date, 6)
    
    # Output
    res <- list(nuc = nuc, loc = loc, date = final.date, nb.snp = length(loc), d.dist = d.dist)
    return(res)
  }
  
  # Sequence evolution by time step
  for (t in time.vector) {
    id <- rownames(cases.info[cases.info$infection.time == t, ])
    snp[id] <- lapply(id, evol, t = t)
  }
  
  ###############################################
  ## WRITING SEQUENCES
  ## Sampling time of sequences
  samp.dates <- sapply(names(snp), function(x) snp[[x]]$date)
  
  # Name of the sequences
  # num.seq_ID_patch_inftime-in-decimal-year
  long.names <- paste(names(snp), samp.dates, sep = "_")
  short.names <- names(snp)
  
  # Location of sequences
  patch <- sapply(long.names, function(x) as.numeric(strsplit(x, "_")[[1]][1])) + 1
  
  # Write text file containing the spatial and temporal metadata
  if (!is.null(coords)) {
    continuous_traits <- data.frame(traits = long.names, lat = coords$lat[patch+1],
                                    long = coords$long[patch+1], regions = coords$regions[patch+1], 
                                    patch = patch, samp.dates = samp.dates, short.names = short.names)
  } else {
    continuous_traits <- data.frame(traits = long.names, patch = patch, 
                                    samp.dates = samp.dates, short.names = short.names)
  }
  
  write.table(continuous_traits, paste0(directory, "/sim", simulationNb, "_metadata.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # Write sequences 
  return(snp)
  
}
