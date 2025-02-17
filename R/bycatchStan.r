#Function to get WAIC and LOOIC from one stan model object
getIC <- function(mod1) {
  LL1 <- extract_log_lik(mod1, "LL")
  waicval <- waic(LL1)$estimates
  looval <- loo(LL1)$estimates
  c(waic = waicval[3, 1], looic = looval[3, 1])
}
#Function to calculate total bycatch from one stan model object
getBycatch <- function(mod1, logdat) {
  gg1 <- extract(mod1, pars = "StrataBycatch")$StrataBycatch
  gg1 <- reshape2::melt(gg1)
  gg1 <- as.data.frame(gg1) %>%
    rename(row = Var2)
  gg1$Year <- logdat$Year[gg1$row]
  modelyrSum1 <- gg1 %>% group_by(Year, iterations) %>%
    summarize(yearsum = sum(value)) %>%
    group_by(Year) %>%
    summarize(
      mean = mean(yearsum),
      se = sd(yearsum),
      lower = quantile(yearsum, 0.025),
      upper = quantile(yearsum, 0.975),
      median = quantile(yearsum, 0.5)
    )
  modelyrSum1
}
#Write out the stan file for negative binomial with estimated effort
write(
  "data{
 int N;
 int Nall;
 int Ncoef;
 int Y[N];
 vector[N] offset;
 matrix[N,Ncoef] xMatrix;
 matrix[Nall,Ncoef] xMatrixAll;
 vector[Nall] EffortMean;
 vector[Nall] EffortSd;
}
parameters{
 vector[Ncoef] b;
 real<lower=0.00001> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  logmu = xMatrix*b;
  for(i in 1:N) {
   mu[i] = exp(logmu[i])*offset[i];
  }
}
model{
  b~normal(0,10);
  phi~normal(0,1);
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  real muAll[Nall];
  real EffortEst[Nall];
  real StrataBycatch[Nall];
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
  for(i in 1:Nall) {
    muAll[i] = exp(xMatrixAll[i,]*b);
    EffortEst[i] =  normal_rng(EffortMean[i],EffortSd[i]);
    StrataBycatch[i] = muAll[i]*EffortEst[i];
  }
}

",
file="NB2matrix.stan")

#Write out the stan file for negative binomial with ordinary logbook effort
write(
  "data{
 int N;
 int Nall;
 int Ncoef;
 int Y[N];
 vector[N] offset;
 matrix[N,Ncoef] xMatrix;
 matrix[Nall,Ncoef] xMatrixAll;
 vector[Nall] EffortAll;
}
parameters{
 vector[Ncoef] b;
 real<lower=0.00001,upper=100> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  logmu = xMatrix*b;
  for(i in 1:N) {
   mu[i] = exp(logmu[i])*offset[i];
  }
}
model{
  b~normal(0,10);
  phi~normal(0,1);
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  real muAll[Nall];
  real StrataBycatch[Nall];
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
  for(i in 1:Nall) {
    muAll[i] = exp(xMatrixAll[i,]*b)*EffortAll[i];
    StrataBycatch[i] = neg_binomial_2_rng(muAll[i],phi);
  }
}


",file="NB2matrixCompleteEffort.stan")

#Function to run a set of negative binomial stan models to estimate bycatch
#taking a bycatchEstimator setup object as an input.
bycatchStan <- function(setupObj,
                        modelsToRun = NULL,
                        spNum = 1,
                        #which of the species to run from multispecies setuObj
                        stanModel = "nbinom2",
                        modeledEffort = FALSE,
                        outDir = NULL) {
  require(rstan)
  require(loo)
  options(mc.cores = parallel::detectCores())
  # To keep a compiled version of the code so you don't have to recompile
  rstan_options(auto_write = TRUE)
  #Unpack setupObj
  modelTry <- obsdat <- logdat <- yearVar <- obsEffort <- logEffort <- logUnsampledEffort <-
    includeObsCatch <- matchColumn <- factorNames <- randomEffects <- randomEffects2 <-
    EstimateIndex <- EstimateBycatch <- logNum <- sampleUnit <- complexModel <-
    simpleModel <- indexModel <-
    designMethods <- designVars <- designPooling <- poolTypes <- pooledVar <-
    adjacentNum <-
    minStrataUnit <-
    baseDir <- runName <- runDescription <-
    common <- sp <- obsCatch <- catchUnit <- catchType <- NULL
  
  numSp <- modelTable <- modelSelectTable <- modFits <- modPredVals <- modIndexVals <-
    residualTab <- bestmod <- predbestmod <- indexbestmod <- allmods <- allindex <-
    modelFail <- rmsetab <- metab <- dat <- yearSum <- requiredVarNames <-
    allVarNames <- indexDat <- strataSum <- NumCores <- NULL
  
  for (r in 1:NROW(setupObj$bycatchInputs))
    assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]])
  for (r in 1:NROW(setupObj$bycatchOutputs))
    assign(names(setupObj$bycatchOutputs)[r], setupObj$bycatchOutputs[[r]])
  
  if (is.null(outDir))
    outDir <- baseDir
  if (!dir.exists(outDir))
    stop("Directory not found")
  numMod <- length(modelsToRun)
  #standardize numeric variables
  numericalVars <- allVarNames[!allVarNames %in% factorNames]
  meanVals <- NULL
  sdVals <- NULL
  if (length(numericalVars) > 0) {
    for (i in 1:length(numericalVars)) {
      meanVals[i] <- mean(obsdat[[numericalVars[i]]], na.rm = TRUE)
      sdVals[i] <- sd(obsdat[[numericalVars[i]]], na.rm = TRUE)
      obsdat[numericalVars[[i]]] <- (obsdat[numericalVars[[i]]] - meanVals[i]) /
        sdVals[i]
      logdat[numericalVars[[i]]] <- (logdat[numericalVars[[i]]] - meanVals[i]) /
        sdVals[i]
    }
  }
  modelTables <- list()
  matrixAll <- list()
  obsdat <- mutate(obsdat, y = 1)
  logdat <- mutate(logdat, y = 1)
  for (i in 1:numMod) {
    mod1 <- lm(formula = formula(modelsToRun[i]), data = obsdat)
    modelTables[[i]] <- model.matrix(mod1)
    matrixAll[[i]] <- model.matrix(formula(mod1), data = logdat)
  }
  stanRunFiles <- NULL
  waicList <- list()
  modelYearSum <- list()
  diagList <- list()
  dirVal <- paste0(outDir, "/Output ", runName, "/", common[spNum], " ", catchType[spNum], "/")
  if (!dir.exists(dirVal))
    dir.create(paste0(dirVal), recursive = TRUE)
  for (i in 1:numMod) {
    if (modeledEffort) {
      dataList <- list(
        Y = obsdat[[obsCatch[spNum]]],
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]),
        offset = obsdat$Effort,
        xMatrix = modelTables[[i]],
        xMatrixAll = matrixAll[[i]],
        EffortMean = logdat$Effort,
        EffortSd = ifelse(logdat$hoursSD > 0, logdat$hoursSD, 0.01 *
                            logdat$Effort),
        Nall = nrow(matrixAll[[i]])
      )
      stanRun <- stan(file = "NB2matrix.stan", data = dataList)
    }  else {
      #If effort is not a model output
      dataList <- list(
        Y = obsdat[[obsCatch[spNum]]],
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]),
        offset = obsdat$Effort,
        xMatrix = modelTables[[i]],
        xMatrixAll = matrixAll[[i]],
        EffortAll = logdat$Effort,
        Nall = nrow(matrixAll[[i]])
      )
      stanRun <- stan(file = "NB2matrixCompleteEffort.stan", data = dataList)
    }
    waicList[[i]] <- getIC(stanRun)
    names(waicList)[i] <- modelsToRun[i]
    modelYearSum[[i]] <- getBycatch(stanRun, logdat = logdat)
    names(modelYearSum)[i] <- modelsToRun[i]
    diagList[[i]] <- data.frame(summary(stanRun, pars = c("b", "phi"))$summary) %>%
      rownames_to_column(var = "Parameter")
    names(diagList)[i] <- modelsToRun[i]
    stanRunFiles[i] <- paste0(dirVal, sp, spNum, "run", i, "-", Sys.Date(), ".rds")
    saveRDS(stanRun, file = stanRunFiles[i])
    rm("stanRun")
  }
  waictab <- bind_rows(waicList, .id = "Model") %>%
    mutate(waic = waic - min(waic), looic = looic - min(looic))
  yearSum <- bind_rows(modelYearSum, .id = "Model")
  diagTable <- bind_rows(diagList, .id = "Model")
  returnVal <- list(
    waictab = waictab,
    yearSum = yearSum,
    diagTable = diagTable,
    stanRunFiles = stanRunFiles
  )
  saveRDS(returnVal, file = paste0(dirVal, Sys.Date(), "StanOutputs.rds"))
  return(returnVal)
}

#Function to plot annual total bycatch from the annual summary table
plotStan <- function(yearSum) {
  ggplot(
    yearSum,
    aes(
      x = Year,
      y = median,
      ymin = lower,
      ymax = upper,
      fill = Model,
      color = Model,
      group = Model
    )
  ) +
    geom_line() +
    geom_ribbon(alpha = 0.2)
}


#write out binomial stan model for survival/mortality. No predictions
write(
  "data{
 int N;
 int Ncoef;
 int Y[N];
 matrix[N,Ncoef] xMatrix;
}
parameters{
 vector[Ncoef] b;
}
transformed parameters{
  vector[N] logitmu;
  logitmu = xMatrix*b;
}
model{
  b~normal(0,10);
  Y~bernoulli_logit(logitmu);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  for(i in 1:N) {
   Yrep[i] = bernoulli_logit_rng(logitmu[i]);
   LL[i] = bernoulli_logit_lpmf(Y[i]|logitmu[i]);
  }
}
", file="R/binomial.stan")

#Write out binomial stan model with new data to predict
write(
  "data{
 int N;
 int Nall;
 int Ncoef;
 int Y[N];
 matrix[N,Ncoef] xMatrix;
 matrix[Nall,Ncoef] xMatrixAll;
}
parameters{
 vector[Ncoef] b;
}
transformed parameters{
  vector[N] logitmu;
  logitmu = xMatrix*b;
}
model{
  b~normal(0,10);
  Y~bernoulli_logit(logitmu);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  real logitMuAll[Nall];
  real strataProb[Nall];
  for(i in 1:N) {
   Yrep[i] = bernoulli_logit_rng(logitmu[i]);
   LL[i] = bernoulli_logit_lpmf(Y[i]|logitmu[i]);
  }
  for(i in 1:Nall) {
    logitMuAll[i] = xMatrixAll[i,]*b;
    strataProb[i] = inv_logit(logitMuAll[i]);
  }
}

", file="R/binomialP.stan")


#Function to run binomial stan models to estimate probability of survival
# 1 is survive, 0 is dead in the alive variable.
mortalityStan <- function(mortData,
                          #Data frame to fit model
                          predData = NULL,
                          #Data frame to predict mortalities for if desired
                          modelsToRun,
                          #Character vector of models
                          aliveColumn,
                          #Name of column containing 1/0 for alive/dead
                          outDir,
                          #output directory
                          runName,
                          #run name
                          predictP) {
  #TRUE/FALSE, do we want to predict to new data?
  
  require(rstan)
  require(loo)
  if (!dir.exists(outDir))
    stop("Directory not found")
  options(mc.cores = parallel::detectCores())
  # To keep a compiled version of the code so you don't have to recompile
  rstan_options(auto_write = TRUE)
  numMod <- length(modelsToRun)
  modelTables <- list()
  matrixAll <- list()
  mortData <- rename(mortData, y = !!aliveColumn)
  predData <- mutate(predData, y = 1)
  for (i in 1:numMod) {
    mod1 <- lm(formula = formula(modelsToRun[i]), data = mortData)
    modelTables[[i]] <- model.matrix(mod1)
    if (predictP)
      matrixAll[[i]] <- model.matrix(formula(mod1), data = predData)
  }
  stanRuns <- list()
  waicList <- list()
  diagList <- list()
  if (predictP)  {
    mod <- stan_model(file = "R/binomialP.stan")
  }  else  {
    mod <- stan_model(file = "R/binomial.stan")
  }
  for (i in 1:numMod) {
    if (predictP) {
      dataList <- list(
        N = nrow(modelTables[[i]]),
        Nall = nrow(matrixAll[[i]]),
        Ncoef = ncol(modelTables[[i]]),
        Y = mortData$y,
        xMatrix = modelTables[[i]],
        xMatrixAll = matrixAll[[i]]
      )
      stanRuns[[i]] <- sampling(mod, data = dataList)
      
    } else  {
      dataList <- list(
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]),
        Y = mortData$y,
        xMatrix = modelTables[[i]]
      )
      stanRuns[[i]] <- sampling(mod, data = dataList)
      
    }
    waicList[[i]] <- getIC(stanRuns[[i]])
    names(waicList)[i] <- modelsToRun[i]
    diagList[[i]] <- data.frame(summary(stanRuns[[i]], pars = c("b"))$summary) %>%
      rownames_to_column(var = "Parameter")
    names(diagList)[i] <- modelsToRun[i]
  }
  waictab <- bind_rows(waicList, .id = "Model") %>%
    mutate(waic = waic - min(waic), looic = looic - min(looic))
  diagTable <- bind_rows(diagList, .id = "Model")
  returnVal <- list(waictab = waictab,
                    diagTable = diagTable,
                    stanRuns = stanRuns)
  dirVal <- paste0(outDir, "/mortality ", runName, "/")
  if (!dir.exists(dirVal))
    dir.create(dirVal, recursive = TRUE)
  saveRDS(returnVal, file = paste0(dirVal, Sys.Date(), "StanOutputsM.rds"))
  return(returnVal)
}

#Function to combine total bycatch with mortality to get total mortality

###

# standardize numeric variables to means and variances from obsdat
# Apply this to logdat so that predictions will be correct if using numerical variables
standardizeToObsdat <- function(obsdat, newdat, numericalVars = NULL) {
  meanVals <- NULL
  sdVals <- NULL
  if (length(numericalVars) > 0) {
    for (i in 1:length(numericalVars)) {
      meanVals[i] <- mean(obsdat[[numericalVars[i]]], na.rm = TRUE)
      sdVals[i] <- sd(obsdat[[numericalVars[i]]], na.rm = TRUE)
      newdat[[paste0("original", numericalVars[i])]] <- newdat[[numericalVars[i]]]
      newdat[[numericalVars[i]]] <- (newdat[[numericalVars[i]]] - meanVals[i]) /
        sdVals[i]
    }
  }
  newdat
}

getMortPred <- function(logdat,
                        mortPredDat,
                        sp,
                        mortMod,
                        bycatchMod) {
  ggb <- extract(bycatchMod, pars = "StrataBycatch")$StrataBycatch %>%
    reshape2::melt() %>%
    as.data.frame() %>%
    rename(row = Var2, StrataBycatch = value) %>%
    mutate(Year = logdat$Year[row])
  ggm <- extract(mortMod, pars = "strataProb")$strataProb %>%
    reshape2::melt() %>%
    as.data.frame() %>%
    rename(row = Var2, strataProb = value) %>%
    mutate(Year = mortPredDat$Year[row], Species = mortPredDat$Species[row]) %>%
    filter(Species == sp) %>%
    mutate(row = row - min(row) + 1)
  gg1 <- left_join(ggb, ggm, by = c("Year", "row", "iterations"))
  modelyrSum1 <- gg1 %>% group_by(Year, iterations) %>%
    summarize(
      Bycatch = sum(StrataBycatch),
      Mortality = sum(StrataBycatch * strataProb)
    ) %>%
    ungroup() %>%
    pivot_longer(Bycatch:Mortality,
                 names_to = "Outcome",
                 values_to = "Total") %>%
    group_by(Year, Outcome) %>%
    summarize(
      mean = mean(Total),
      se = sd(Total),
      lower = quantile(Total, 0.025),
      upper = quantile(Total, 0.975),
      median = quantile(Total, 0.5)
    )
  modelyrSum1
}

#Function to plot both bycatch and bycatch mortality
plotMortalityFunc <- function(modelyrSum1, Species) {
  ggplot(
    modelyrSum1,
    aes(
      x = Year,
      y = median,
      fill = Outcome,
      color = Outcome,
      ymin = lower,
      ymax = upper
    )
  ) +
    geom_line() +
    geom_ribbon(alpha = 0.4) +
    scale_fill_manual(values = c("blue", "red")) +
    scale_color_manual(values = c("blue", "red")) +
    labs(
      y = paste0("Total number of ", Species),
      fill = "",
      color = ""
    )
}




#Write out the stan file for negative binomial with
#no simulation of effort
write(
  "data{
 int N;
 int Ncoef;
 int Y[N];
 vector[N] offset;
 matrix[N,Ncoef] xMatrix;
}
parameters{
 real  b0;
 vector[Ncoef] b;
 real<lower=0.00001,upper=100> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  logmu = b0+xMatrix*b;
  for(i in 1:N) {
   mu[i] = exp(logmu[i])*offset[i];
  }
}
model{
  b0~normal(0,10);
  b~normal(0,1);
  phi~normal(0,1);
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
}

",file="NB2matrixNoEffort.stan")

#No effort intercept only
write(
  "data{
 int N;
 int Ncoef;
 int Y[N];
 vector[N] offset;
}
parameters{
 real  b0;
 real<lower=0.00001,upper=100> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  for(i in 1:N) {
   logmu[i] = b0;
   mu[i] = exp(logmu[i])*offset[i];
  }
}
model{
  b0~normal(0,10);
  phi~normal(0,1);
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
}

",file="NB2matrixNoEffort1.stan")


#Function to run a set of negative binomial stan models to estimate bycatch
#taking a bycatchEstimator setup object as an input. and using simulation
bycatchStanSim <- function(setupObj,
                           modelsToRun = NULL,
                           spNum = 1,
                           #which of the species to run from multispecies setuObj
                           stanModel = "nbinom2",
                           modeledEffort = FALSE,
                           effortSD = NULL,
                           predictionInterval=TRUE,
                           outDir = NULL) {
  if (is.null(effortSD) & modeledEffort)      stop("Must supply the name of the effortSD column if using estimated effort")
  require(rstan)
  require(loo)
  options(mc.cores = parallel::detectCores())
  # To keep a compiled version of the code so you don't have to recompile
  rstan_options(auto_write = TRUE)
  #Unpack setupObj
  modelTry <- obsdat <- logdat <- yearVar <- obsEffort <- logEffort <- logUnsampledEffort <-
    includeObsCatch <- matchColumn <- factorNames <- randomEffects <- randomEffects2 <-
    EstimateIndex <- EstimateBycatch <- logNum <- sampleUnit <- complexModel <-
    simpleModel <- indexModel <-
    designMethods <- designVars <- designPooling <- poolTypes <- pooledVar <-
    adjacentNum <-
    minStrataUnit <-
    baseDir <- runName <- runDescription <-
    common <- sp <- obsCatch <- catchUnit <- catchType <- NULL
  
  numSp <- modelTable <- modelSelectTable <- modFits <- modPredVals <- modIndexVals <-
    residualTab <- bestmod <- predbestmod <- indexbestmod <- allmods <- allindex <-
    modelFail <- rmsetab <- metab <- dat <- yearSum <- requiredVarNames <-
    allVarNames <- indexDat <- strataSum <- NumCores <- NULL
  
  for (r in 1:NROW(setupObj$bycatchInputs))
    assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]])
  for (r in 1:NROW(setupObj$bycatchOutputs))
    assign(names(setupObj$bycatchOutputs)[r], setupObj$bycatchOutputs[[r]])
  
  if (is.null(outDir))
    outDir <- baseDir
  if (!dir.exists(outDir))
    stop("Directory not found")
  numMod <- length(modelsToRun)
  #standardize numeric variables
  numericalVars <- allVarNames[!allVarNames %in% factorNames]
  meanVals <- NULL
  sdVals <- NULL
  obsdat$Year1 <- obsdat$Year
  logdat$Year1 <- logdat$Year
  if (length(numericalVars) > 0) {
    for (i in 1:length(numericalVars)) {
      if(is.numeric(meanVals[i])) {
       meanVals[i] <- mean(obsdat[[numericalVars[i]]], na.rm = TRUE)
       sdVals[i] <- sd(obsdat[[numericalVars[i]]], na.rm = TRUE)
       obsdat[numericalVars[[i]]] <- (obsdat[numericalVars[[i]]] - meanVals[i]) /
         sdVals[i]
       logdat[numericalVars[[i]]] <- (logdat[numericalVars[[i]]] - meanVals[i]) /
         sdVals[i]
      }
     }
  }
  modelTables <- list()
  matrixAll <- list()
  obsdat <- mutate(obsdat, y = 1)
  logdat <- mutate(logdat, y = 1)
  for (i in 1:numMod) {
    mod1 <- lm(formula = formula(modelsToRun[i]), data = obsdat)
    modelTables[[i]] <- model.matrix(mod1)
    matrixAll[[i]] <- model.matrix(formula(mod1), data = logdat)
  }
  stanRunFiles <- NULL
  waicList <- list()
  modelYearSum <- list()
  diagList <- list()
  dirVal <- paste0(outDir, "/Output ", runName, "/", common[spNum], " ", catchType[spNum], "/")
  if (!dir.exists(dirVal))
    dir.create(paste0(dirVal), recursive = TRUE)
  for (i in 1:numMod) {
    if (modelsToRun[i] == "y~1") {
      dataList <- list(
        Y = obsdat[[obsCatch[spNum]]],
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]) - 1,
        offset = obsdat$Effort
      )
      stanRun <- stan(file = "NB2matrixNoEffort1.stan", data = dataList,
                      pars=c("b0","phi","LL"))
    } else {
      dataList <- list(
        Y = obsdat[[obsCatch[spNum]]],
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]) - 1,
        offset = obsdat$Effort,
        xMatrix = as.matrix(modelTables[[i]][, -1])
      )
      stanRun <- stan(file = "NB2matrixNoEffort.stan", data = dataList,
                      pars=c("b0","b","phi","LL"))
    }
    waicList[[i]] <- getIC(stanRun)
    names(waicList)[i] <- modelsToRun[i]
    modelYearSum[[i]] <- getBycatchSim(
      stanRun,
      logdat = logdat,
      matrixAll = matrixAll[[i]],
      modeledEffort = modeledEffort,
      effortSD = effortSD,
      predictionInterval=predictionInterval,
      
    )
    names(modelYearSum)[i] <- modelsToRun[i]
    if (modelsToRun[i] == "y~1")
      pars <- c("b0", "phi")
    else
      pars <- c("b0", "b", "phi")
    diagList[[i]] <- data.frame(summary(stanRun, pars = pars)$summary) %>%
      rownames_to_column(var = "Parameter")
    names(diagList)[i] <- modelsToRun[i]
    stanRunFiles[i] <- paste0(dirVal, sp, spNum, "run", i, "-", Sys.Date(), ".rds")
    saveRDS(stanRun, file = stanRunFiles[i])
    rm("stanRun")
  }
  waictab <- bind_rows(waicList, .id = "Model") %>%
    mutate(waic = waic - min(waic), looic = looic - min(looic))
  yearSum <- bind_rows(modelYearSum, .id = "Model")
  diagTable <- bind_rows(diagList, .id = "Model")
  returnVal <- list(
    waictab = waictab,
    yearSum = yearSum,
    diagTable = diagTable,
    stanRunFiles = stanRunFiles
  )
  saveRDS(returnVal, file = paste0(dirVal, Sys.Date(), "StanOutputs.rds"))
  return(returnVal)
}

#Function for a random draw of size SampleUnits, summed to get
#the stratum estimate. 
getMeanNbinom<-function(SampleUnits,MeanVals,phiVals) {
  if(SampleUnits>0 & !is.na(SampleUnits) & !is.na(MeanVals))
    return<- sum(rnbinom(SampleUnits,mu=MeanVals/SampleUnits,size=phiVals)) else
      return<-0
    return
}
#Function to calculate total bycatch from one stan model object
#simulating the catches in R not stan, with prediction interval generated with GetMeanNbinom
getBycatchSim <- function(mod1,
                          logdat,
                          matrixAll,
                          modeledEffort = FALSE,
                          effortSD = NULL,
                          predictionInterval=predictionInterval,
                          nsim = 1000,
                          usePrior=FALSE,
                          returnDraws=FALSE) {
  if(usePrior) b0vals<-rnorm(nsim,0,10) else
   b0vals <- extract(mod1, pars = "b0")$b0
  subsetval <- sample(1:length(b0vals), nsim)
  if(ncol(matrixAll) > 1)
    bvals <- extract(mod1, pars = "b")$b else
      bvals <- NULL
  if(!is.null(bvals) & usePrior)
     bvals[,]<-rnorm(prod(dim(bvals)))
  b0vals <- b0vals[subsetval]
  bvals <- bvals[subsetval, ]
  bvals <- cbind(b0vals, bvals)
  phivals <- extract(mod1, pars = "phi")$phi[subsetval]
  simMean <- exp(matrixAll %*% t(bvals))
  EffortMean <- rep(logdat$Effort, nsim)
  EffortMean[EffortMean < 0.01] <- 0.01
  if(modeledEffort) {
    EffortSD <- rep(unlist(logdat[, effortSD]), nsim)
    Effort <- rnorm(n = prod(dim(simMean)), EffortMean, EffortSD)  #Check this
    Effort[Effort < 0.0001] <- 0.0001
  } else {
    Effort <- EffortMean
  }
  if(predictionInterval) {
  simVal<-data.frame(SampleUnits=rep(logdat$SampleUnits,nsim),
                     MeanVals=as.vector(simMean) * Effort,
                     phiVals=rep(phivals, each = nrow(logdat))) %>%
    rowwise() %>%
    mutate(Bycatch=getMeanNbinom(SampleUnits,MeanVals,phiVals))
  }  else {
    simVal<-data.frame(Bycatch=as.vector(simMean) * Effort)
  }
  # simVal <- rnbinom(
  #   n = prod(dim(simMean)),
  #   mu = as.vector(simMean) * Effort,
  #   size = rep(phivals, each = nrow(logdat))
  # )
  gg1 <- data.frame(
    simVal = simVal$Bycatch,
    row = rep(1:nrow(logdat),nsim),
    iterations = rep(1:nsim, each = nrow(logdat))
  )
  gg1$Year <- logdat$Year1[gg1$row]
  modelyrSum1 <- gg1 %>% group_by(Year, iterations) %>%
    summarize(yearsum = sum(simVal)) %>%
    group_by(Year) %>%
    summarize(
      mean = mean(yearsum),
      se = sd(yearsum),
      lower = quantile(yearsum, 0.025),
      upper = quantile(yearsum, 0.975),
      median = quantile(yearsum, 0.5)
    )
  return<-modelyrSum1
  if(returnDraws) return<-list(yearSum=modelyrSum1,
                               gg1=gg1)
}

# prior simulation from default priors
priorSimulation<-function(stanObj,coefs,nsim=1000) {
  return<-data.frame(Iteration=1:nsim,
                     Chain=1,
                     Parameter=rep(coefs,each=nsim)) %>%
    rowwise() %>%
    mutate(value=case_when(Parameter=="b0"~rnorm(1,0,10),
                           Parameter=="phi"~rexp(1,1),
                           TRUE~rnorm(1,0,1)))
  return
}

# plot prior and posterior
plotPriorPosterior<-function(stanObj) {
  posterior<-ggs(stanObj) %>%
    filter(grepl("b",Parameter) | Parameter=="phi") 
  coefs<-as.character(unique(posterior$Parameter))
  prior<-posteriorSimulation(stanObj,coefs)
  df<-bind_rows(list(prior=prior,posterior=posterior),.id="type")
  ggplot(df,aes(x=value,color=type))+
    geom_density()+
    facet_wrap(~Parameter,scales="free")  +
    scale_color_manual(values=c("blue","darkgrey"))+
    labs(x="Parameter",y="Density",color="")
}
