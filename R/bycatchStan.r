library(tidyverse)
library(BycatchEstimator)
library(MuMIn)
library(gridExtra)
library(rstan)
library(loo)
library(shinystan)
library(ggmcmc)
library(readxl)
library(bayesplot)
library(ggsci)
library(flextable)
theme_set(theme_bw())

#' getIC
#' Function to extract WAIC and LOOIC from fitted stan model object
#'
#'
#' @param mod1 rstan or cmndstanr object
#' @param useCode value "cmdstanr" or "rstan" to indicate code used in fit
#'
getIC <- function(mod1,useCode) {
  if(useCode=="rstan")  LL1 <- extract_log_lik(mod1, "LL")
  if(useCode=="cmdstanr") LL1 <-mod1$draws("LL")
  waicval <- waic(LL1)$estimates
  looval <- loo(LL1)$estimates
  c(waic = waicval[3, 1], looic = looval[3, 1])
}

#' plotStan
#' Function to plot annual total bycatch from the annual summary table
#'
#' @param yearSum annual bycatch summary table 
#'
#' @export
#' 
plotStan <- function(yearSum) {
  print(ggplot(
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
  )
}


#' mortalityStan
#' Function to run binomial stan models to estimate probability of survival
#'
#' @param mortData Data frame to fit model with a variable called alive with 1 for survive, 0 for dead plus predictor varialbes.
#' @param predData Data frame to predict mortalities for if desired
#' @param modelsToRun  Character vector of models
#' @param aliveColumn Name of column containing 1/0 for alive/dead
#' @param outDir output directory
#' @param runName run name
#' @param predictP TRUE/FALSE, do we want to predict to new data?
#' @param useCode character string: cmdstanr or rstan
#'
#' @export
#' 
mortalityStan <- function(mortData,
                          predData = NULL,
                          modelsToRun,
                          aliveColumn,
                          outDir,
                          runName,
                          predictP,  
                          useCode  
                          ) {
  if(useCode=="rstan") {
    require(rstan)
    rstan_options(auto_write = TRUE)
  }  
  if(useCode=="cmdstanr") {
    require(cmdstanr)
    binomialMod<-cmdstan_model("stan/binomial.stan")
    binomialModP<-cmdstan_model("stan/binomialP.stan")
  }
  require(loo)
  if (!dir.exists(outDir))
    stop("Directory not found")
  options(mc.cores = parallel::detectCores())
  numMod <- length(modelsToRun)
  modelTables <- list()
  matrixAll <- list()
  mortData <- rename(mortData, y = !!aliveColumn)
  if(!is.null(predData))predData <- mutate(predData, y = 1)
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
    if(useCode=="rstan") mod <- stan_model(file = "R/binomialP.stan")
    if(useCode=="cmdstanr") mod <- binomialModP
  }  else  {
    if(useCode=="rstan") mod <- stan_model(file = "R/binomial.stan")
    if(useCode=="cmdstanr") mod <- binomialMod  
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
    } else  {
      dataList <- list(
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]),
        Y = mortData$y,
        xMatrix = modelTables[[i]]
      )
    }
    if(useCode=="rstan") stanRuns[[i]] <- sampling(mod, data = dataList)
    if(useCode=="cmdstanr") stanRuns[[i]] <-mod$sample(data=dataList,
                                                       refresh=0)
    waicList[[i]] <- getIC(stanRuns[[i]],useCode)
    names(waicList)[i] <- modelsToRun[i]
    if(useCode=="rstan") diagList[[i]] <- data.frame(summary(stanRuns[[i]], pars = c("b"))$summary) %>%
      rownames_to_column(var = "Parameter")
    if(useCode=="cmdstanr") diagList[[i]] <- stanRuns[[i]]$summary(variables = "b") %>%
      dplyr::rename(Parameter = variable)
    names(diagList)[i] <- modelsToRun[i]
  }
  waictab <- bind_rows(waicList, .id = "Model") %>%
    mutate(Dwaic = waic - min(waic), Dlooic = looic - min(looic))
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

#' standardizeToObsdat
#' Function to standardize numeric variables to means and variances from obsdat
#' Apply this to logdat so that predictions will be correct if using numerical variables
#'
#' @param obsdat observer data used in model fitting
#' @param newdat data to be standardized, e.g. logbook data
#' @param numericVariables character vector of names of numeric variables
#'
#' @export
#' 
standardizeToObsdat <- function(obsdat, newdat, numericVariables = NULL) {
  meanVals <- NULL
  sdVals <- NULL
  if (length(numericVariables) > 0) {
    for (i in 1:length(numericVariables)) {
      meanVals[i] <- mean(obsdat[[numericVariables[i]]], na.rm = TRUE)
      sdVals[i] <- sd(obsdat[[numericVariables[i]]], na.rm = TRUE)
      newdat[[paste0("original", numericVariables[i])]] <- newdat[[numericVariables[i]]]
      newdat[[numericVariables[i]]] <- (newdat[[numericVariables[i]]] - meanVals[i]) /
        sdVals[i]
    }
  }
  newdat
}

#' plotMortalityFunc
#' Function to plot both bycatch and bycatch mortality
#' @param modelyrSum1 Dataframe with data to plot
#' @param Species Common name of species
#'
#' @returns Prints out a ggplot
#' @export
#'
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


#' bycatchStanSim
#' Function to run a set of negative binomial stan models to estimate bycatch
#' taking a bycatchEstimator setup object as an input, and using simulation for effort if needed
#'
#' @param setupObj List output from a rund of BycatchEstimator::bycatchSetup
#' @param modelsToRun Character vector of models to run, e.g. c("y~Year","y~1")
#' @param spNum Number of the species in the original bycatchSetup run, 
#' @param stanModel Type of likelihood to use, currently only negative binomial 2
#' @param priors List of priors. 
#' @param modeledEffort TRUE/FALSE for whether effort is normally distributed with a mean and SD
#' @param effortSD  Name of effortSD variable in logdat if used. 
#' @param predictionInterval TRUE to esimate prediction interval rather than confidence interval
#' @param StanOutDir Directory for output. NULL defaults to same directory as BycatchEstimator outputs
#' @param useCode "stanr" or "cmdstanr"
#'
#' @returns Returns lists of all the inputs, as well as a model summary table with LOOIC and WAIC values
#' and a vector with paths of the .rds files containing the individual model stan files, and a dataframe of 
#' annual bycatch estimates suitable for plotting
#' @export
#'
#' @examples
bycatchStanSim <- function(setupObj,
                           modelsToRun = NULL,
                           spNum = 1,
                           #which of the species to run from multispecies setuObj
                           stanModel = "nbinom2",
                           priors =  list(interceptSD=10,
                                          coefficientSD=1,
                                          phiType=c("exponential","normal")[1],
                                          phiPar=1),
                           modeledEffort = FALSE,
                           effortSD = NULL,
                           predictionInterval=TRUE,
                           StanOutDir = NULL,
                           useCode="cmdstanr") {
  if (is.null(effortSD) & modeledEffort)      stop("Must supply the name of the effortSD column if using estimated effort")
  if(useCode=="rstan") require(rstan)
  if(useCode=="cmdstanr") {
    require(cmdstanr)
    NB2matrixNoEffort <- cmdstan_model("stan/NB2matrixNoEffort.stan")
    NB2matrixNoEffort1 <- cmdstan_model("stan/NB2matrixNoEffort1.stan")
  }
  if(!useCode %in% c("rstan","cmdstanr")) stop("Must specify rstan or cmdstanr")
  require(loo)
  options(mc.cores = parallel::detectCores())
  # To keep a compiled version of the code so you don't have to recompile
  if(useCode=="rstan") rstan_options(auto_write = TRUE)
  #Unpack setupObj
  obsdat<-logdat<-yearVar<-obsEffort<-logEffort<-obsCatch<-catchUnit<-catchType<-
    logNum<-sampleUnit<-factorVariables<-numericVariables<-EstimateBycatch<-
    baseDir<-dirname<-outDir<-runName<-runDescription<-common<-sp<-NULL
   dat<-numSp<-yearSum<-allVarNames<-startYear<-strataSum<-shortName<-NULL
   for(r in 1:NROW(setupObj$bycatchInputs)) assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]])
   if(all(is.na(numericVariables))) numericVariables<-NULL
  
  if(is.null(StanOutDir))
    StanOutDir <- outDir
  if(!dir.exists(StanOutDir))
    stop("Directory not found")
  numMod <- length(modelsToRun)
  #standardize numeric variables
  meanVals <- NULL
  sdVals <- NULL
  obsdat$Year1 <- obsdat$Year
  logdat$Year1 <- logdat$Year
  if(length(numericVariables) > 0) {
    for (i in 1:length(numericVariables)) {
      if(is.numeric(obsdat[,meanVals[i]])) {
       meanVals[i] <- mean(obsdat[[numericVariables[i]]], na.rm = TRUE)
       sdVals[i] <- sd(obsdat[[numericVariables[i]]], na.rm = TRUE)
       obsdat[numericVariables[[i]]] <- (obsdat[numericVariables[[i]]] - meanVals[i]) /
         sdVals[i]
       logdat[numericVariables[[i]]] <- (logdat[numericVariables[[i]]] - meanVals[i]) /
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
  dirVal<-paste0(StanOutDir, "/",abbreviate(common[spNum],minlength = 10),
                                 abbreviate(catchType[spNum],minlength = 5),"/STAN")
  if(!dir.exists(dirVal))
    dir.create(paste0(dirVal), recursive = TRUE)
  for (i in 1:numMod) {
    if(modelsToRun[i] == "y~1") {
      dataList <- list(
        Y = obsdat[[obsCatch[spNum]]],
        N = nrow(modelTables[[i]]),
        interceptSD=priors$interceptSD,
        phiType=ifelse(priors$phiType=="normal",2,1),
        phiPar=priors$phiPar,
        Ncoef = ncol(modelTables[[i]]) - 1,
        Effort = obsdat$Effort
      )
      if(useCode=="rstan") {
        stanRun <- stan(file = "stan/NB2matrixNoEffort1.stan", data = dataList,
                        pars=c("b0","phi","LL"))
      }
      if(useCode=="cmdstanr")  {
        stanRun<-NB2matrixNoEffort1$sample(data=dataList,
                                           refresh=0)
      }
    } else {
      dataList <- list(
        Y = obsdat[[obsCatch[spNum]]],
        N = nrow(modelTables[[i]]),
        Ncoef = ncol(modelTables[[i]]) - 1,
        Effort = obsdat$Effort,
        xMatrix = as.matrix(modelTables[[i]][, -1]),
        interceptSD = priors$interceptSD,
        coefficientSD = priors$coefficientSD,
        phiType=ifelse(priors$phiType=="normal",2,1),
        phiPar=priors$phiPar
      )
      if(useCode=="rstan") {
        stanRun <- stan(file = "stan/NB2matrixNoEffort.stan", data = dataList,
                        pars=c("b0","b","phi","LL"))
      }
      if(useCode=="cmdstanr")  {
        stanRun<-NB2matrixNoEffort$sample(data=dataList,
                                          refresh=0)
      }
    }
    waicList[[i]] <- getIC(stanRun,useCode)
    names(waicList)[i] <- modelsToRun[i]
    modelYearSum[[i]] <- getBycatchSim(
      stanRun,
      logdat = logdat,
      matrixAll = matrixAll[[i]],
      modeledEffort = modeledEffort,
      effortSD = effortSD,
      priors =  priors,
      predictionInterval=predictionInterval,
      useCode=useCode
    )
    names(modelYearSum)[i] <- modelsToRun[i]
    if (modelsToRun[i] == "y~1")
      pars <- c("b0", "phi") else
      pars <- c("b0", "b", "phi")
    if(useCode=="rstan") 
      diagList[[i]] <- data.frame(summary(stanRun, pars = pars)$summary) %>%
        rownames_to_column(var = "Parameter")
    if(useCode=="cmdstanr") diagList[[i]] <- stanRun$summary(variables = pars)
    names(diagList)[i] <- modelsToRun[i]
    stanRunFiles[i] <- paste0(dirVal,"/", sp, spNum, "run", i, "-", Sys.Date(), ".rds")
    #print(stanRunFiles[i])
    if(useCode=="rstan") saveRDS(stanRun, file = stanRunFiles[i])
    if(useCode=="cmdstanr") stanRun$save_object(file = stanRunFiles[i])
    rm("stanRun")
  }
  waictab <- bind_rows(waicList, .id = "Model") %>%
    mutate(Dwaic = waic - min(waic), Dlooic = looic - min(looic))
  yearSum <- bind_rows(modelYearSum, .id = "Model")
  diagTable <- bind_rows(diagList, .id = "Model")
  returnVal <- list(
    waictab = waictab,
    yearSum = yearSum,
    diagTable = diagTable,
    stanRunFiles = stanRunFiles,
    stanInputs = list(modelsToRun=modelsToRun,
                      spNum=spNum,
                      modeledEffort=modeledEffort,
                      effortSD=effortSD,
                      predictionInterval=predictionInterval,
                      priors=priors)
  )
  saveRDS(returnVal, file = paste0(dirVal, Sys.Date(), "StanOutputs.rds"))
  return(returnVal)
}

#' getMeanNbinom
#' Function for a random draw of size SampleUnits, summed to get the stratum estimate of bycatch. 
#'
#' @param SampleUnits 
#' @param MeanVals 
#' @param phiVals 
#'
getMeanNbinom<-function(SampleUnits,MeanVals,phiVals) {
  if(SampleUnits>0 & !is.na(SampleUnits) & !is.na(MeanVals))
    return<- sum(rnbinom(SampleUnits,mu=MeanVals/SampleUnits,size=phiVals)) else
      return<-0
    return
}

#' getBycatchSim
#' Function to calculate total bycatch from one stan model object
#' simulating the catches in R not stan, with prediction interval generated with GetMeanNbinom
#'
#' @param mod1 stan fit object from cmdstanr or rstan
#' @param logdat Total effort data for expanding over
#' @param matrixAll Model matrix from total effort data
#' @param modeledEffort TRUE for effort being drawn from a normal with mean and SD, false for effort input by sample unit
#' @param effortSD Optional SD for effort if drawing froma  noraml
#' @param predictionInterval TRUE for prediction interval, false for confidence interval
#' @param nsim Number of draws for calculation
#' @param usePrior TRUE to simulate from priors instead of using fit
#' @param priors Priors used in the model fit
#' @param returnDraws TRUE to return random draws
#' @param useCode "cmdstanr" or "rstan"
#'
#' @returns  Returns lists of inputs and outputs
#' @export
#'
#' @examples
getBycatchSim <- function(mod1,
                          logdat,
                          matrixAll,
                          modeledEffort = FALSE,
                          effortSD = NULL,
                          predictionInterval=predictionInterval,
                          nsim = 1000,
                          usePrior=FALSE,
                          priors = list(interceptSD=10,
                                        coefficientSD=1,
                                        phiType=c("exponential","normal")[1],
                                        phiPar=1),
                          returnDraws=FALSE,
                          useCode="cmdstanr") {
  if(usePrior) {
    b0vals<-rnorm(nsim,0,priors$interceptSD) 
    if(ncol(matrixAll)>1)
      bvals<-matrix(rnorm(prod(nsim,(ncol(matrixAll)-1)),0,priors$coefficientSD),nsim,(ncol(matrixAll)-1))
    if(priors$phiType=="normal") 
      phivals=truncnorm::rtruncnorm(nsim, a=0, b=Inf, mean = 0, sd = priors$phiPar) else
    phivals<-rexp(nsim,priors$phiPar)
  }
  else {
    if(useCode=="rstan") b0vals <- extract(mod1, pars = "b0")$b0
    if(useCode=="cmdstanr") b0vals <- mod1$draws(variables="b0",format = "df")$b0
    subsetval <- sample(1:length(b0vals), nsim)
    if(ncol(matrixAll) > 1) {
      if(useCode=="rstan") bvals <- extract(mod1, pars = "b")$b
      if(useCode=="cmdstanr") bvals <- mod1$draws(variables="b",format="matrix")
    } else
        bvals <- NULL
    b0vals <- b0vals[subsetval]
    bvals <- bvals[subsetval, ]
    if(useCode=="rstan") phivals <- extract(mod1, pars = "phi")$phi[subsetval]
    if(useCode=="cmdstanr") phivals <- mod1$draws(variables="phi",format="df")$phi[subsetval] 
  }
  bvals <- cbind(b0vals, bvals)
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
  return
}

#' priorSimulation
#' Prior simulation from specified priors
#'
#' @param stanObj Fitted cmdstanr or rstan object
#' @param coefs Names of coefficients to simulate
#' @param nsim Number of simulations
#' @param priors List of priors
#'
priorSimulation<-function(stanObj,
                          coefs,
                          nsim=1000,
                          priors =  list(interceptSD=10,
                                         coefficientSD=1,
                                         phiType=c("exponential","normal")[1],
                                         phiPar=1)
) {
  return<-data.frame(Iteration=1:nsim,
                     Chain=1,
                     Parameter=rep(coefs,each=nsim)) %>%
    rowwise() %>%
    mutate(value=case_when(Parameter=="b0"~rnorm(1,0,priors$interceptSD),
                           Parameter=="phi"&priors$phiType=="exponential"~rexp(1,priors$phiPar),
                           Parameter=="phi"&priors$phiType=="normal"~truncnorm::rtruncnorm(1,0,Inf,0,priors$phiPar),
                           TRUE~rnorm(1,0,priors$coefficientSD)))
  return
}

# plot prior and posterior
#' plotPriorPosterior
#'
#' @param stanSum Output of bycatchStanSim
#' @param stanObj Fitted cmdstanr or rstan object
#' @param useCode "cmdstanr" or "rstan"
#'
#' @returns
#' @export
#'
#' @examples
plotPriorPosterior<-function(stanSum, 
                             stanObj,
                             useCode) {
  if(useCode=="rstan") posterior<-ggs(stanObj) %>%
    filter(grepl("b",Parameter) | Parameter=="phi") 
  if(useCode=="cmdstanr") posterior<-ggs(as_mcmc.list(stanObj))%>%
      filter(grepl("b",Parameter) | Parameter=="phi")
  coefs<-as.character(unique(posterior$Parameter))
  prior<-priorSimulation(stanObj,
                         coefs,
                         nsim=1000,
                         priors = stanSum$stanInputs$priors )
  df<-bind_rows(list(prior=prior,posterior=posterior),.id="type")
  ggplot(df,aes(x=value,color=type))+
    geom_density()+
    facet_wrap(~Parameter,scales="free")  +
    scale_color_manual(values=c("blue","darkgrey"))+
    labs(x="Parameter",y="Density",color="")
}

#' plotPriorPosteriorSims
#' A function to plot Prior and Posterior predictive checks
#'
#' @param stanSum Object output by `bycatchStanSim`
#' @param modelNum Which model number to look at (matches waictab rows)
#' @param stanObj The rstan or cmdstanr fit object
#' @param setupObj The BycatchEstimator setup object
#' @param modeledEffort TRE if effort is drawn from a normal
#' @param effortSD Effort SD column in obsdat if drawing from noraml
#' @param useCode "cmdstanr" or "rstan"
#'
#' @returns
#' @export
#'
#' @examples
plotPriorPosteriorSims<-function(stanSum,
                                 modelNum,
                                 stanObj,
                                 setupObj,
                                 modeledEffort = FALSE,
                                 effortSD = NULL,
                                 useCode)  {
  logdat<-setupObj$bycatchInput$logdat %>% 
    mutate(y=1,
           Year1=Year)
  matrixAll<-model.matrix(formula(stanSum$waictab$Model[modelNum]),data=logdat)
  postvals<-getBycatchSim(stanObj,
                          logdat,
                          matrixAll,
                          modeledEffort = modeledEffort,
                          effortSD = effortSD,
                          predictionInterval=TRUE,
                          nsim = 30,
                          usePrior=FALSE,
                          returnDraws=TRUE,
                          useCode=useCode)[[2]] %>% 
    group_by(Year, iterations) %>%
    summarize(yearsum = sum(simVal)) 
  
  priorvals<-getBycatchSim(stanFit,
                           logdat,
                           matrixAll,
                           modeledEffort = modeledEffort,
                           effortSD = effortSD,
                           predictionInterval=TRUE,
                           nsim = 30,
                           usePrior=TRUE,
                           returnDraws=TRUE,
                           useCode=NULL)[[2]] %>% 
    group_by(Year, iterations) %>%
    summarize(yearsum = sum(simVal)) 
  g1<-ggplot(priorvals,aes(x=as.numeric(as.character(Year)),y=yearsum))+
    geom_line(aes(group=iterations),alpha=0.9,color="lightblue")+
    geom_line(data=filter(stanSum$yearSum,Model==stanSum$waictab$Model[modelNum]),
              aes(x=as.numeric(as.character(Year)),y=mean),color="darkred")+
    labs(x="Year",y="Bycatch",color="",title="Annual bycatch prior draws")
  g2<-ggplot(postvals,aes(x=as.numeric(as.character(Year)),y=yearsum))+
    geom_line(aes(group=iterations),alpha=0.9,color="lightblue")+
    geom_line(data=filter(stanSum$yearSum,Model==stanSum$waictab$Model[modelNum]),
              aes(x=as.numeric(as.character(Year)),y=mean),color="darkred")+
    labs(x="Year",y="Bycatch",color="",title="Annual bycatch posterior draws")
  gridExtra::grid.arrange(g1,g2)
}

#' getResiduals
#' Function to plot quantile residuals from a cmdstanr or rstan fit, using DHARMa
#'
#' @param stanSum Output from bycatchStanSim
#' @param stanFit One rstan or cmdstanr object
#' @param modelNum Model number corresponding to rows in waictab
#' @param setupObj bycatchSetup output
#' @param spNum  Species number from bycatchSetup, generally 1
#' @param useCode "cmdstanr" or "rstan"

#' @returns
#' @export
#'
#' @examples
getResiduals<-function(stanSum,stanFit,modelNum,setupObj,nsim=1000,spNum=1,useCode) {
  #require(rstan)
  require(DHARMa)
  obsdat<-setupObj$bycatchInputs$obsdat
  obsdat$y<-1
  obsdat$y<-obsdat[[setupObj$bycatchInputs$obsCatch[spNum]]]
  if(useCode=="rstan") b0vals <- extract(stanFit, pars = "b0")$b0
  if(useCode=="cmdstanr") b0vals <- stanFit$draws(variables="b0",format="df")$b0
  subsetval <- sample(1:length(b0vals), nsim)
  matrix1<-model.matrix(formula(stanSum$waictab$Model[modelNum]),data=obsdat)
  if(ncol(matrix1) > 1) {
      if(useCode=="rstan") bvals <- extract(stanFit, pars = "b")$b
      if(useCode=="cmdstanr") bvals <- stanFit$draws(variables="b",format="matrix")
      bvals <- bvals[subsetval, ]
  } else
      bvals <- NULL
  b0vals <- b0vals[subsetval]
  bvals <- cbind(b0vals, bvals)
  if(useCode=="rstan") phivals <- extract(stanFit, pars = "phi")$phi[subsetval]
  if(useCode=="cmdstanr") phivals <- stanFit$draws(variables="phi",format="df")$phi[subsetval] 
  simMean <- exp(matrix1 %*% t(bvals))
  EffortMean <- rep(obsdat$Effort, nsim)
  EffortMean[EffortMean < 0.01] <- 0.01
  simVal <- rnbinom(
    n = prod(dim(simMean)),
    mu = as.vector(simMean) * EffortMean,
    size = rep(phivals, each = nrow(obsdat))
  )
  simVal<-matrix(simVal,nrow(obsdat),nsim)
  DHARMaRes <- createDHARMa(simulatedResponse =simVal , 
                            observedResponse = obsdat$y, 
                            fittedPredictedResponse = apply(simVal,1,mean),
                            integerResponse = TRUE)
  plot(DHARMaRes)
}

#' getSummary
#' Function to extract parameter summary table for coefficients and scale parameter
#'
#' @param stanSum Output from bycatchStanSim
#' @param stanFit One rstan or cmdstanr object
#' @param modelNum Model number corresponding to rows in waictab
#' @param setupObj bycatchSetup output
#' @param spNum  Species number from bycatchSetup, generally 1
#' @param useCode "cmdstanr" or "rstan"
#'
#' @returns
#' @export
#'
#' @examples
getSummary<-function(stanSum,stanFit,modelNum,setupObj,spNum=1,useCode) {
  obsdat<-setupObj$bycatchInputs$obsdat
  obsdat$y<-1
  obsdat$y<-obsdat[[setupObj$bycatchInputs$obsCatch[spNum]]]
  matrixAll<-model.matrix(formula(stanSum$waictab$Model[modelNum]),data=obsdat)
  parVals<-c("phi","b0")
  if(ncol(matrixAll)>1)
    parVals<-c(parVals,paste0("b[",1:(ncol(matrixAll)-1),"]"))
  if(useCode=="rstan") temp<-summary(stanObj,pars=parVals)$summary
  if(useCode=="cmdstanr") temp<-stanObj$summary(variables=parVals)
  temp
}

#' Get convergence diagnostics
#'
#' @param stanObj Fitted cmdstanr or rstan object
#' @param useCode "cmdstanr" or "rstan"
#'
#' @returns Rstan convergence plots or cmdstanr number of diverenges,etc.
#' @export
#'
#' @examples
getConvergence<-function(stanObj,useCode) {
 if(useCode=="rstan") print(stan_diag(stanObj))
 if(useCode=="cmdstanr") {
   temp<-stanObj$diagnostic_summary()
   print(temp)
 }
}


#Function to get mortality predictions 
#' getMortPred
#'
#' @param logdat 
#' @param mortPredDat 
#' @param sp 
#' @param mortMod 
#' @param bycatchMod 
#' @param codeName 
#'
#' @returns
#' @export
#'
#' @examples
getMortPred <- function(logdat,
                        mortPredDat,
                        sp,
                        mortMod,
                        bycatchMod,
                        codeName) {
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

