# install.packages("devtools")
# devtools::install_github("ebabcock/BycatchEstimator", force=TRUE)
library(tidyverse)
library(BycatchEstimator)
library(MuMIn)
library(gridExtra)
library(rstan)
library(cmdstanr)
library(loo)
library(shinystan)
library(ggmcmc)
library(readxl)
library(bayesplot)
library(ggsci)
library(flextable)
theme_set(theme_bw())

source("R/bycatchStan.r")

## Run example code from bycatchEstimator

setupObjBUM<-bycatchSetup(
  obsdat = droplevels(LLSIM_BUM_Example_observer[LLSIM_BUM_Example_observer$Year>2010 &LLSIM_BUM_Example_observer$fleet==2,]),
  logdat = droplevels(LLSIM_BUM_Example_logbook[LLSIM_BUM_Example_logbook$Year>2010 & LLSIM_BUM_Example_logbook$fleet==2,]),
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  factorVariables = c("Year","area"),  
  numericVariables = NA, 
  EstimateBycatch = TRUE,
  logNum = NA,
  sampleUnit = "trips",
  baseDir = getwd(),
  runName = "LLSIMBUMtripExample",
  runDescription = "LLSIm BUM by trip, with 5% observer coverage",
  common = "Blue marlin",
  sp = "Makaira nigricans",
  obsCatch ="BUM",
  catchUnit = "number",
  catchType = "catch",
  reportType = "html"  
)
modelsToRun<-c("y~Year","y~1")
BUMRun<-bycatchStanSim(setupObjBUM,
                        modelsToRun=modelsToRun,
                        spNum=1,  #which species to run in input is multispecies
                        modeledEffort=FALSE,
                        priors =  list(interceptSD=10,
                                      coefficientSD=1,
                                      phiType=c("exponential","normal")[1],
                                      phiPar=1),
                        outDir=getwd(),
                        useCode=c("cmdstanr","rstan")[1])
Sys.time()
#Function prints to an rds file. You can read it in here.

#Plot bycatch estimates
plotStan(BUMRun$yearSum)

#Look at model selection table
BUMRun$waictab

# Code to reload runs to look at diagnostics
#Specify directory with bycatch Estimator results
outDir<-paste0(getwd(),"/Output LLSIMBUMtripExample/")
# Date of bycatch estimator run
estimatorDate<-"2025-07-02"
# Date of stan run
stanDate<-"2025-07-02"
# Read in bycatchEstimator setup
setupObj<-readRDS(c(paste0(outDir,estimatorDate,"_BycatchSetupSpecification.rds")))
spNum<-1
# Read in stan run summary files
stanSum<-readRDS(paste0(outDir,setupObj$bycatchInputs$common[spNum]," ",setupObj$bycatchInputs$catchType,
                        "/",stanDate,"StanOutputs.rds"))
stanSum$waictab
modelNum<-which(stanSum$waictab$waic==0) #To get WAIC best
# Load the selected model
stanObj<-readRDS(paste0(outDir,setupObj$bycatchInputs$common[spNum]," ",setupObj$bycatchInputs$catchType,
                        "/",setupObj$bycatchInputs$sp[spNum],spNum,"run",modelNum,"-",
                        stanDate,".rds"))

# Plot residuals
getResiduals(stanSum, 
             stanObj,
             modelNum,
             setupObj)

# plot prior and posterior densities of parameters
plotPriorPosterior(stanObj) 

# Plot prior and posterior draws of bycatch 
plotPriorPosteriorSims(stanSum,
                       modelNum,
                       stanObj,
                       setupObj,
                       modeledEffort = FALSE,
                       effortSD = NULL)  

#Convergence diagnostics
stan_diag(stanObj)
getSummary(stanSum,stanFit,modelNum,setupObj,spNum=1) 

