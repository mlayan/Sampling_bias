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
surveillance <- c('uniformS_10', 'maxPerRegion_10', 'maxPerRegionYear_10',
                  'uniformS_20', 'maxPerRegion_20', 'maxPerRegionYear_20')

bias2 <- c('uniform', 'biased_10', 'biased_20')
surveillance2 <- c('uniformS', 'maxPerRegionYear')

# All sampling protocols
all = c('uniform_150', 'biased_2.5_150', 'biased_5_150', 'biased_10_150', 'biased_20_150', 'biased_50_150', 
        'uniform_500', 'biased_2.5_500', 'biased_5_500', 'biased_10_500', 'biased_20_500', 'biased_50_500',
        'uniformS_150', 'maxPerRegion_150', 'maxPerRegionYear_150', 
        'uniformS_500', 'maxPerRegion_500', 'maxPerRegionYear_500')

# Regions
# regions = c("WesternMediterraneancoast", "EasternMediterraneancoast", 
#             "NorthAtlanticcoast", "SouthAtlanticcoast", "NorthernAtlas", 
#             'SouthernAtlas', "AntiAtlas")

regions = c("Region1", "Region2", "Region3", "Region4", "Region5", 'Region6', "Region7")

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
  s <- as.numeric(param[2])
  p <- param[3]
  if (length(param) == 4) ma <- param[4]
  
  # Subset the dataframe
  if (exists("ma")) {
    D = inputData[inputData$nSim == si &
                    inputData$nSeq == s &
                    inputData$protocol == p & 
                    inputData$matrix == ma, ]
  } else {
    D = inputData[inputData$nSim == si & 
                    inputData$nSeq == s &
                    inputData$protocol == p, ]
  }
  
  # Location frequencies
  f <- D[D$parameter == "locationFrequency" & order(D$source), c('value', 'source')]
  
  # Diagonal of the Q matrix
  ratesDiag <- D %>%
    filter(source != dest, parameter == "migrationRate") %>% 
    group_by(source) %>% 
    summarise(qdiag = sum(value, na.rm = TRUE), .groups = 'drop') %>% 
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
  p <- as.character(param[2]) # Protocol
  if (length(param) == 3) {  
    # Matrix 
    ma <- param[3]
    if (ma == "1" | ma == "2") ma = as.numeric(ma)
  } 
  
  print(paste(s,p))
  
  # Median colname
  if (stat == "median") {
    if ("X50" %in% colnames(df)) medianCol = "X50"
    if ("median" %in% colnames(df)) medianCol = "median"
  }
  
  # BEAST models tested
  models = as.character(unique(df$model[df$nSeq == s & df$protocol == p]))
  
  # Output
  out = data.frame()
  
  for (model in models) {
    
    # Subset dataframe 
    df_sub = df[df$nSeq == s & df$protocol == p & df$model == model, ]
    if (length(param) == 3) {
      df_sub <- df_sub[df_sub$matrix == ma, ] 
    }
    
    if (nrow(df_sub) == 0 ) {
      if (test == "lm") outTemp = data.frame(nSeq = s, 
                                             protocol = p,  
                                             model = model,
                                             coef = NA, 
                                             r = NA, 
                                             slope = NA, 
                                             intercept = NA)
      
      if (test == "kendall") outTemp = data.frame(nSeq = s, 
                                                  protocol = p,  
                                                  model = model,
                                                  coef = NA, 
                                                  tau = NA, 
                                                  p = NA)
      if (test == "spearman") outTemp = data.frame(nSeq = s, 
                                                   protocol = p,  
                                                   model = model,
                                                   coef = NA, 
                                                   rho = NA, 
                                                   p = NA) 
      
      out = bind_rows(out, outTemp)
      
    } else {
      
      # Compute regression coefficient
      if (test == "lm") {
        if (stat == "median") lModel <- lm(df_sub[[medianCol]] ~ df_sub$value.y)
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
                              coef = paste0(as.character(as.expression(eq))), 
                              r = coef, 
                              slope = slope, 
                              intercept = intercept) 
        
        out = bind_rows(out, outTemp)
        
      } else if (test == "spearman") {
        if (stat == 'mean') sModel <- cor.test(df_sub$mean, df_sub$value.y, method = "spearman",alternative ="two.sided")
        if (stat == 'median') sModel <- cor.test(df_sub[[medianCol]], df_sub$value.y, method = "spearman",alternative ="two.sided")
        pVal = sModel$p.value
        r = sModel$estimate
        
        eqS <- paste0("rho == ~\'", format(r, digits = 3),"\' - p ==~\'",format(pVal, digits = 3), "\'") 
        # paste(as.expression(
        # substitute(rho==~rh-p==~pV, 
        #            list(rh=format(r, digits = 3), pV=format(pVal, digits = 3))
        # )))
        #paste0("rho '=' ", format(sModel$estimate, digits = 3), 
        #      " - p '=' ", format(sModel$p.value, digits = 3))
        
        
        # Append results to output
        outTemp <- data.frame(nSeq = s, 
                              protocol = p,  
                              model = model,
                              coef = eqS, 
                              rho = r, 
                              p = pVal) 
        out = bind_rows(out, outTemp)
        
      } else if (test == "kendall") {
        if (stat == 'mean') sModel <- cor.test(df_sub$mean, df_sub$value.y, method = "kendall", alternative ="two.sided")
        if (stat == 'median') sModel <- cor.test(df_sub[[medianCol]], df_sub$value.y, method = "kendall", alternative ="two.sided")
        pVal = sModel$p.value
        t = sModel$estimate
        
        eqK <- paste0("tau == ~\'", format(t, digits = 3),"\' - p ==~\'",format(pVal, digits = 3), "\'") 
        #paste0("tau '=' ", format(sModel$estimate, digits = 3), 
        #      " - p '=' ", format(sModel$p.value, digits = 3))
        
        # Append results to output
        outTemp <- data.frame(nSeq = s, 
                              protocol = p,  
                              model = model,
                              coef = eqK, 
                              tau = t, 
                              p = pVal) 
        out = bind_rows(out, outTemp)
        
      } else {
        warning("The test argument is not valid")
      }
    }
  }
  
  if (length(param) == 3) out$matrix = ma
  return(out)
  
}
#############################################################
# Weighted interval scores
WIS = function(data) {
  
  require(scoringutils)
  
  K = length(data)
  
  IS = sapply(data, function(x) {
    interval_score(
      true_values = x[["smooth_value"]],
      lower=x[["lower"]],
      upper=x[["upper"]],
      interval_range=unique(x[["ir"]]),
      weigh = TRUE,
      separate_results = FALSE)
  })
  
  IS_100<-abs(data[[1]]$point-data[[1]]$smooth_value)*0.5
  wis = (apply(IS, 1, sum)+IS_100)/(K+0.5)
  return(wis)
}
