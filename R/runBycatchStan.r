# install.packages("devtools")
# devtools::install_github("ebabcock/BycatchEstimator")
library(tidyverse)
library(BycatchEstimator)
library(MuMIn)
library(gridExtra)
library(rstan)
library(loo)
library(shinystan)
library(ggmcmc)
library(readxl)
library(lubridate)
library(bayesplot)
library(ggsci)
library(flextable)

source("R/bycatchStan.r")

## Run example code from bycatchEstimator

setupObjBUM<-bycatchSetup(
  modelTry = c("TMBnbinom1"),
  obsdat = droplevels(LLSIM_BUM_Example_observer[LLSIM_BUM_Example_observer$Year>2010 &LLSIM_BUM_Example_observer$fleet==2,]),
  logdat = droplevels(LLSIM_BUM_Example_logbook[LLSIM_BUM_Example_logbook$Year>2010 & LLSIM_BUM_Example_logbook$fleet==2,]),
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  logUnsampledEffort = "unsampledEffort",
  factorNames = c("Year","area"),
  EstimateBycatch = TRUE,
  logNum = NA,
  sampleUnit = "trips",
  complexModel = formula(y~Year+area),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year),
  baseDir = getwd(),
  runName = "LLSIMBUMtripExample",
  runDescription = "LLSIm BUM by trip, with 5% observer coverage",
  common = "Blue marlin",
  sp = "Makaira nigricans",
  obsCatch ="BUM",
  catchUnit = "number",
  catchType = "catch"
)

dataCheck(setupObjBUM)
modelsToRun<-c("y~Year+area","y~Year","y~1")
BUMRun<-bycatchStanSim(setupObjBUM,
                        modelsToRun=modelsToRun,
                        spNum=1,  #which species to run in input is multispeces
                        modeledEffort=FALSE,
                        outDir=getwd())
Sys.time()
#Function prints to an rds file. You can read it in here.
#BUMRUn<-readRDS("Output LLSIMBUMtripExample/Blue Marlin Catch/2025-01-06StanOutputs.rds")

#Plot bycatch estimates
plotStan(BUMRun$yearSum)

#Look at model selection table
BUMRun$waictab

# Read in stan run for best model and look at diagnostics
bestMod<-readRDS("Output LLSIMBUMtripExample/Blue Marlin Catch/Makaira nigricans1run1-2025-01-06.rds")
stan_dens(bestMod,par="b")
stan_diag(bestMod)

