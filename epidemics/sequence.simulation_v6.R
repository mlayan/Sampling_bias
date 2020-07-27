#############################################################
##                  SIMULATE SEQUENCES  
##                        
##
## Creation Date: 19-07-2019
## Last Update: 07-10-2019
## Maylis Layan
#############################################################


sequence.simulation <- function(sim, ref, mut.rate, root.date, simulationNb, coords = NULL, bases.comp = NULL, genome.length = NULL) {
  
  #################################################
  ## PACKAGES AND PARAMETERS
  # Load packages
  require(seqinr, quietly = TRUE)
  require(lubridate, quietly = TRUE)
  
  # Directory
  directory <- paste0("simulation", simulationNb, "/")
  
  # Vector of bases
  bases <- c("A","C","G","T")
  
  # Convert cases.info from a list to a data frame
  cases.info <- data.frame(do.call("rbind", sim$sourcing$data))
  names(cases.info) <- c("patch.source", "ID.source", "infection.time", "death.time", "infectious.time", 
                         "infectious.period", "generation.time", "dead", "hard.stop")
  
  # Number of initial sequences
  n.init <- length(which(cases.info$infection.time < 0))
  names.init <- rownames(cases.info[which(cases.info$infection.time < 0),])
  
  # Generation of the sequence of reference
  # By default, it is the first sequence appearing in cases.info
  if (!is.null(ref)) {
    genome.length <- length(ref)
    bases.comp <- as.numeric(table(ref))
  }
  
  if (is.null(bases.comp)) bases.comp <- c(0.25,0.25,0.25,0.25)
  if (is.null(genome.length)) genome.length <- 11937
  if (is.null(ref)) ref <- sample(bases, genome.length, prob = bases.comp, replace = TRUE)
  
  # Eliminate individuals who die before becoming infectious
  cases.info <- cases.info[cases.info$infectious.time >= 0, ]
  
  # Lists of SNPs 
  snp <- vector("list", nrow(cases.info))
  names(snp) <- rownames(cases.info)
  snp[[names.init[1]]] <- list(loc = c(), nuc = ref, date = decimal_date(as.Date(root.date)), nb.snp = 0, d.dist = 0)
  
  ################################################
  ## INDEX CASES
  if (n.init > 1) {
    for (i in 2:n.init) {
      loc <- sample(1:genome.length, round(0.1*genome.length, 0))
      nuc <- sapply(loc, function(x) sample(bases[!(bases %in% ref[x])],1))
      snp[[names.init[i]]] <- list(loc = loc, nuc = nuc, date = root.date, nb.snp = 0, d.dist = 0)
    }
  }
  
  # Time max can exceed the time.max used for the simulation if the epidemic
  # doesn't stop before the end of the simulation 
  # To speed up the evolution loop, we define only the time steps that 
  # should be treated
  time.vector <- sort(unique(cases.info$infection.time))
  time.vector <- time.vector[which(time.vector >= 0)] 
  
  ################################################
  ## SEQUENCE EVOLUTION
  # Evolution function
  evol <- function(x, t) {
    vec <- cases.info[rownames(cases.info) == x, ]
    generation.time <- as.numeric(vec[7])
    
    nb.snp <- rpois(1, mut.rate * genome.length * generation.time)
    loc <- sample(1:genome.length, nb.snp)
    
    nuc <- c()
    ID.infector <- paste(vec[1:2], collapse = "_")
    
    if (ID.infector == names.init[1]) {
      
      nuc <- sapply(loc, function(x) sample(bases[!(bases %in% ref[x])], 1))
      
    } else {
      
      
      nuc <- sapply(loc, function(x) { ifelse(x %in% snp[[ID.infector]]$loc,
                                              sample(bases[!(bases %in% snp[[ID.infector]]$nuc[x==snp[[ID.infector]]$loc])], 1),
                                              sample(bases[!(bases %in% ref[x])], 1))
      })
      
      # Add the SNPs of the infector
      new.snp <- which(snp[[ID.infector]]$loc %in% loc)
      if (length(new.snp)) {
        nuc <- c(nuc, snp[[ID.infector]]$nuc[-new.snp])
        loc <- c(loc, snp[[ID.infector]]$loc[-new.snp]) 
      } else {
        nuc <- c(nuc, snp[[ID.infector]]$nuc)
        loc <- c(loc, snp[[ID.infector]]$loc) 
      }
      
    }
    
    # root.date should be expressed as "yyyy-mm-dd", ie "%Y-%m-%d" or "%Y/%m/%d"
    final.date <- decimal_date(as.Date(root.date) + as.numeric(vec[5]))
    
    # Output
    res <- list(nuc = nuc, loc = loc, date = final.date, nb.snp = nb.snp, d.dist = length(nuc))
    return(res)
  }
  
  # Sequence evolution by time step
  for (t in time.vector) {
    
    id <- rownames(cases.info[cases.info$infection.time == t, ])
    snp[id] <- lapply(id, function(x) evol(x,t))
    
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
  
  write.table(continuous_traits, paste0(directory, "sim", simulationNb, "_metadata.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # Write sequences 
  return(snp)
  
}

