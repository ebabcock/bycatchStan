# install.packages("devtools")
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

source("R/printStanCode.r")
source("R/bycatchStan.r")

## Run example code from BycatchEstimator

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
useCode="rstan"
BUMRun<-bycatchStanSim(setupObjBUM,
                        modelsToRun=modelsToRun,
                        spNum=1,  #which species to run in input is multispecies
                        modeledEffort=FALSE,
                        priors =  list(interceptSD=10,
                                      coefficientSD=1,
                                      phiType=c("exponential","normal")[1],
                                      phiPar=1),
                        StanOutDir=NULL,  #Null to use base directory from SetupObj
                        useCode=c("cmdstanr","rstan")[2])
Sys.time()
#Function prints to an rds file. You can read it in here.

#Plot bycatch estimates
plotStan(BUMRun$yearSum)

#Look at model selection table
BUMRun$waictab

##########################################################
# Code to reload runs to look at diagnostics
# Not needed if the objects are still in the environment
#################################################
#Specify directory with bycatch Estimator results
#outDir<-setupObjBUM$bycatchInputs$outDir
# Date of bycatch estimator run
#estimatorDate<-Sys.Date()
# # Date of stan run
# stanDate<-Sys.Date()
# # Read in bycatchEstimator setup
# setupObj<-readRDS(paste0(outDir,"/",estimatorDate,"_BycatchSetupSpecification.rds"))
# spNum<-1
# #Stan output directory
# StanOutDir<-paste0(outDir,"/","Bluemarlincatch")
# # Read in stan run summary files
# stanSum<-readRDS(paste0(StanOutDir,
#                         "/",stanDate,"StanOutputs.rds"))
#########################################################################

####Using run that was just run
stanSum<-BUMRun
setupObj<-setupObjBUM
############

## Look at WAIC and LOOIC
stanSum$waictab
## Pick AIC best model
modelNum<-which(stanSum$waictab$Dwaic==0) #To get WAIC best
## Read in the Stan model object (these are deleted from the environment to save space)
stanObj<-readRDS(BUMRun$stanRunFiles[[modelNum]])

## Plot residuals
getResiduals(stanSum, 
             stanObj,
             modelNum,
             setupObj,
             useCode=useCode)

# plot prior and posterior densities of parameters
plotPriorPosterior(stanSum, 
                   stanObj,
                   useCode=useCode) 

# Plot prior and posterior draws of bycatch 
plotPriorPosteriorSims(stanSum,
                       modelNum,
                       stanObj,
                       setupObj,
                       modeledEffort = FALSE,
                       effortSD = NULL,
                       useCode=useCode)  

#Convergence diagnostics
getConvergence(stanObj,
               useCode=useCode)

#Parameter summary table
getSummary(stanSum,
           stanFit,
           modelNum,
           setupObj,
           spNum=1,
           useCode=useCode) 


