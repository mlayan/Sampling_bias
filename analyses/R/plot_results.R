#############################################################
##              FUNCTIONS AND DATA TO ANALYZE 
##                  BEAST1 AND BEAST2 RUNS
##
## Creation Date: 19-05-2020
## Last Update: -05-2020
## Maylis Layan
#############################################################


#############################################################
## Simulation parameters
#############################################################
## Clock rate used to simulate sequences
clockRate <- 2.44e-4

## Root age used to simulate sequences
root.age.initial <- 1989

## Kappa 
kappa = 2

## Base frequencies
freqA = 0.2851752
freqC = 0.2197719
freqG = 0.2312594
freqT = 0.2637934 

# ESS threshold
threshold = 100

# Sampling protocols
bias <- c('uniform', 'biased_2.5', 'biased_5', 'biased_10', 'biased_20', 'biased_50')
surveillance <- c('uniformS', 'maxPerRegion', 'maxPerRegionYear')

# All sampling protocols
all = c('uniform_150', 'biased_2.5_150', 'biased_5_150', 'biased_10_150', 'biased_20_150', 'biased_50_150', 
        'uniform_500', 'biased_2.5_500', 'biased_5_500', 'biased_10_500', 'biased_20_500', 'biased_50_500',
        'uniformS_150', 'maxPerRegion_150', 'maxPerRegionYear_150', 
        'uniformS_500', 'maxPerRegion_500', 'maxPerRegionYear_500')

# Regions
regions = c("WesternMediterraneancoast", "EasternMediterraneancoast", 
            "NorthAtlanticcoast", "SouthAtlanticcoast", "NorthernAtlas", 
            'SouthernAtlas', "AntiAtlas")

#############################################################
##  FUNCTIONS
#############################################################
# Compute RMSE
rmse <- function(xhat, x) {
  
  out <- (xhat - x)^2
  out <- sqrt(mean(out))
  
  return(out)
}


#############################################################
RMSE <- function(estimations, paramValues, param) {
  
  # Fixed parameters
  d <- data.frame(pValues = c(freqA, freqC, freqG, freqT, kappa, clockRate*10e4), 
                  param = c('freqA', 'freqC', 'freqG', 'freqT', 'kappa', 'clockRate'))
  
  # Vector of parameters
  p <- unique(param)
  if (length(p) != 1) {
    stop("Only one parameter can be given at a time")
  }
  
  # Compute RMSE
  if (p %in% c('freqA', 'freqC', 'freqG', 'freqT', 'kappa', 'clockRate')) {
    if (p == "clockRate") {
      estimations = estimations * 10e4
    }
    
    out <-rmse(estimations, d$pValues[d$param == p])
    
  } else {
    
    out <- rmse(estimations, paramValues)
    
  }
  
  return(out)
}


#############################################################
# Normalize transition rates 
normalizeRates <- function(param, inputData) {
  require(dplyr)
  
  # Retrieve parameters
  si <- as.numeric(param[1])
  s <- param[2]
  p <- param[3]
  ma <- param[4]
  
  # Subset the dataframe
  if (is.na(s)) {
    D = inputData[inputData$nSim == si & 
                    inputData$protocol == p & 
                    inputData$matrix == ma, ]
  } else {
    D = inputData[inputData$nSim == si & 
                    inputData$nSeq == as.numeric(s) &
                    inputData$protocol == p & 
                    inputData$matrix == ma, ]
  }
  
  # Location frequencies
  f <- D[D$parameter == "locationFrequency" & order(D$source), c('value', 'source')]
  
  # Diagonal of the Q matrix
  ratesDiag <- D %>%
    filter(source != dest, parameter == "migrationRate") %>% 
    group_by(source) %>% 
    summarise(qdiag = sum(value, na.rm = TRUE)) %>% 
    right_join(., f, by = c("source" = "source")) %>%
    data.frame()
  ratesDiag$qdiag[is.na(ratesDiag$qdiag)] <- 0
  
  # Normalizing constant
  normCons <- sum(ratesDiag$qdiag*ratesDiag$value)
  
  out <- D[D$parameter == "migrationRate", ]
  out$value <- out$value / normCons
  return(out)
  
}


#############################################################
# Regression coefficient
test_eqn <- function(param, df, stat = 'median', test = "lm"){
  
  # Retrieve parameters
  s <- as.numeric(param[1])  # Sample size
  p <- param[2]             # Protocol
  if (length(param) == 3) {  
    # Matrix 
    ma <- param[3]
    if (ma == "1" | ma == "2") ma = as.numeric(ma)
  } 
  
  # BEAST models tested
  models = unique(df$model)
  
  # Output
  out = data.frame()
  
  for (model in models) {
    # Subset dataframe 
    df_sub = df[df$nSeq == s & df$protocol == p & df$model == model, ]
    if (length(param) == 3) {
      df_sub <- df_sub[df_sub$matrix == ma, ] 
    }
    
    # Compute regression coefficient
    if (test == "lm") {
      if (stat == "median") lModel <- lm(X50 ~ value.y, df_sub)
      if (stat == "mean") lModel <- lm(mean ~ value.y, df_sub)
      eq <- substitute(#y == a + b %.% x*","~~
        r^2~"="~r2, 
        list(#a = format(unname(coef(m_dta)[1]), digits = 2),
          #b = format(unname(coef(m_dta)[2]), digits = 2),
          r2 = format(summary(lModel)$r.squared, digits = 3)))
      coef = summary(lModel)$r.squared
      intercept = coef(lModel)[1]
      slope = coef(lModel)[2] 
      
      # Append results to output
      outTemp <- data.frame(nSeq = s, 
                            protocol = p,  
                            model = model,
                            coef = as.character(as.expression(eq)), 
                            r = coef, 
                            slope = slope, 
                            intercept = intercept) 
      
      out = bind_rows(out, outTemp)
      
    } else if (test == "spearman") {
      if (stat == 'mean') sModel <- cor.test(df_sub$mean, df_sub$value.y, method = "spearman")
      if (stat == 'median') sModel <- cor.test(df_sub$X50, df_sub$value.y, method = "spearman")
      eqS <- paste0(model, ": ", rho, ' = ', format(sModel$estimate, digits = 3), 
                         ", p = ", format(sModel$p.value, digits = 3))
      pVal = sModel$p.value
      rho = sModel$estimate
      
      # Append results to output
        out <- data.frame(nSeq = s, 
                          protocol = p,  
                          model = model,
                          coef = eqS, 
                          rho = rho, 
                          p = pVal) 
        
    
    } else {
      warning("The test argument is not valide")
    }
  }
  
  if (length(param) == 3) out$matrix = ma
  return(out)
  
}

#############################################################
# Regression coefficient
test_direction <- function(param, df, test = "lm"){
  
  # Retrieve parameters
  s <- as.numeric(param[1])
  p <- param[2]
  m <- param[3]

  # Subset dataframe
  forw <- df[df$nSeq == s & df$protocol == p & df$matrix == m & df$direction == "forwards", ] 
  backw <- df[df$nSeq == s & df$protocol == p & df$matrix == m & df$direction == "backwards", ]

  # Compute regression coefficient
  if (test == "lm") {
    
    m_f <- lm(mean ~ value.y, forw)
    eq_f <- substitute(r^2~"="~r2, 
      list(r2 = format(summary(m_f)$r.squared, digits = 3)))
    coefF = summary(m_f)$r.squared
    interceptF = coef(m_f)[1]
    slopeF = coef(m_f)[2] 
    
    # Compute regression coefficient
    m_b <- lm(mean ~ value.y, backw)
    eq_b <- substitute(r^2~"="~r2, 
                            list(r2 = format(summary(m_b)$r.squared, digits = 3)))
    coefB = summary(m_b)$r.squared
    interceptB = coef(m_b)[1]
    slopeB = coef(m_b)[2] 
    
    # Output as a dataframe
    out <- data.frame(nSeq = rep(s,2), protocol = rep(p, 2), 
                      direction = c("backwards", 'forwards'),
                      coef = c(paste0("backwards: ", as.character(as.expression(eq_b))),
                               paste0("forwards: ", as.character(as.expression(eq_f)))), 
                      r = c(coefB, coefF), 
                      slope = c(slopeB, slopeF), 
                      intercept = c(interceptB, interceptF)) 
    
    
  } else if (test == "spearman") {
    
    s_f <- cor.test(forw$mean, forw$value.y, method = "spearman")
    eq_s_f <- paste0('forwards: = ', format(s_f$estimate, digits = 3), 
                       ", p = ", format(s_f$p.value, digits = 3))
    pF = s_f$p.value
    rhoF = s_f$estimate
    
    # Compute regression coefficient
    s_b <- cor.test(backw$mean, backw$value.y, method = "spearman")
    eq_s_b <- paste0('mascot: r = ', format(s_b$estimate, digits = 3), 
                          ", p = ", format(s_b$p.value, digits = 3))
    pB = s_b$p.value
    rhoB = s_b$estimate
    
    # Output as a dataframe
    out <- data.frame(nSeq = rep(s,2), protocol = rep(p, 2), 
                      direction = c("forwards", 'backwards'),
                      coef = c(eq_s_f, eq_s_b), 
                      rho = c(rhoF, rhoB), 
                      p = c(pF, pB)) 
    
  } else {
    warning("The test argument is not valid")
  }
  
  out$matrix = m
  
  return(out)
}


#############################################################
# Subset dataframe
drop_rows <- function(param, df = data) {
  m = as.numeric(param[1])
  n = as.numeric(param[2])
  s = as.numeric(param[3])
  p = as.character(param[4])
  md = as.character(param[5])
  
  out = subset(df, !(matrix == m && nSim == n && matrix == s && protocol == p && model == md))
  
  return(out)
}


#############################################################
# Paired two-sided wilcoxon test
tests <- function(x) {
  require(dplyr)
  
  seq = as.numeric(x[1])
  mat = as.numeric(x[2])
  prot = as.character(x[3])
  
  d = filter(topologies_wide, protocol == prot, nSeq == seq, matrix == mat)
  out = wilcox.test(d$dta, d$mascot, paired = TRUE, alternative = "two.sided")
  
  return(c(seq, mat, prot, out$p.value))
}

