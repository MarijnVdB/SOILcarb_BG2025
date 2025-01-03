# This is the main script to run the SOILcarb model, as described in
# Van de Broek et al., A microbially-driven and depth-explicit soil organic 
# carbon model constrained by carbon isotopes to reduce equifinality,
# Biogeosciences

# These codes are published under a Creative Commons Attribution-NonCommercial-ShareAlike
# (CC BY-NC-SA) License (https://creativecommons.org/licenses/by-nc-sa/4.0/)

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")

startTime <- proc.time()

# Packages are loaded
library("readxl")
library("ggplot2")
library("deSolve")
library("psych")
library("gridExtra")
library("GA")
library("parallel")
library("doParallel")
library("foreach")
library("iterators")
library("DEoptim")
library("rgenoud")
library("ggpubr")

# The files containing the functions and other scripts are sourced
source("f_various_functions.R")
source("f_load_measured_data.R")
source("f_inputParameters.R")
source("f_load_d13C_data.R")
source("f_load_d14C_data.R")
source("f_construct_isotope_inputs.R")
source("f_initialCO2DepthProfile.R")
source("f_parameter_initialisation.R")
source("f_runCarbonSimulations.R")
source("f_SOILcarb.R")
source("f_plotting.R")
source("f_createClasses.R")
source("AutoCluster4.R")

# ------------------------------------------------------------------------------
# General model parameters are defined 
# ------------------------------------------------------------------------------

# An S4 class object to store the parameters in is created
createInputClass() # Create the function to create the input object
createVegClass()   # Create the function to create the vegetation object
par <- new("inputObject") # An object to store all parameter values in

# The name of the current site
par@site <- "Hainich"

# The number of vegetation stages, and how long they last
par@nVeg = 1      # The number of vegetation stages simulated
par@yearsStage1 <- 15000  # Number of years for vegetation species 1
par@yearsStage2 <- 1000   # Number of years for vegetation species 2: not used here

# The number of total simulation years
if(par@nVeg == 1){
  par@numberOfSimulationYears <- par@yearsStage1
} else if (par@nVeg == 2){
  par@numberOfSimulationYears <- par@yearsStage1 + par@yearsStage2
}

# The last 'real' year of the simulations (calendar year)
if(par@site == "Hainich"){
  par@lastRealSimulationYear <- 2004
}

# Parameters to set the depth layers
layer1_thickness <- 0.01     # Thickness of the uppermost layer [m]
maxDepth_layers <- 1          # The approximate maximum soil depth to be simulated [m]
boundaryFactor <- 0.5

# The layer boundaries and mid depths are calculated
out <- calculateLayers(layer1_thickness, maxDepth_layers, boundaryFactor)
par@boundaryDepth <- out$boundaryDepth
par@midDepth_layer <- out$midDepth
par@layerThickness <- out$layerThickness
rm(out, boundaryFactor, layer1_thickness, maxDepth_layers)

par@totalDepth<- max(par@boundaryDepth)         # The depth over which the profile is simulated [m]
par@layerNumber <- length(par@midDepth_layer)   # The number of simulated layers

# Processes that need to be simulated
par@Include_Suess_effect <- 1               # If 1, the 'terrestrial Suess effect' is simulated
par@includeBioturbation <- 1                # If 1, bioturbation is simulated
par@includeAdvection <- 1                   # If 1, advection is simulated

# ------------------------------------------------------------------------------
# Calibration and sensitivity analysis options
# ------------------------------------------------------------------------------

par@normalRun <-  1                       # If 1, the mode is run in 'normal mode'
par@calibMode <- 0                        # If 1, a frequentist calibration is performed
par@useCalibrationResults <- 1            # If 1, the calibration results are used as input parameters
par@runSeries <- 0                        # To run the model for a series of parameter values (used to evaluate the calibration results)
par@counter_loadIsotopes <- 0             # Necessary in calibration mode (do not change)
par@SAFE_sensitivity_isotopes <- 0        # The parameter combinations for a SAFE sensitivity analysis are run
par@SAFE_sensitivity_allParameters <- 0   # The parameter combinations for a SAFE sensitivity analysis are run

par@evalDepth <- 0.5                      # The depth down to which the C profile is evaluated during calibration

par@parallelCalib <- 1                    # If 1, a parallel calibration is performed

par@litterCalibration <- 0                # If 1, the litter layer is calibrated
par@soilCalibration <- 1                  # If 1, the soil is calibrated

# Some potential errors are identified
if(par@litterCalibration == 1 && par@soilCalibration == 1){
  stop("Error: you cannot run a litter and soil calibration simultaneously!")
}

if(par@calibMode == 1 & par@normalRun == 1){
  stop("Error: you cannot run a calibration and a normal run simultaneously!")
}

if(par@calibMode == 1 & par@useCalibrationResults == 1){
  stop("Error: you cannot calibrate and use the calibrated data as the same time!")
}

# Decide which variables are calibrated during the DE calibration
if(par@calibMode == 1){
  # A parameter to indicate which variables are used for calibration
  par@calibCheck <- data.frame(
    "C" = 0,              # Scenario with C only
    "C_d13C" = 0,         # Scenario with C and d13C
    "C_d14C" = 0,         # Scenario with C and d14C
    "C_d13C_d14C" = 1)    # Scenario with C, d13C and d14C
  if(sum(par@calibCheck > 1)){
    stop("Error: variables that have to be calibrated during DE calibration are not correctly defined!")
  }
}

# Use the calibration results as input
if(par@useCalibrationResults == 1){
  
  # filename = "DE output/C/DE_out.RData"
  # filename = "DE output/C_d13C/DE_out.RData"
  # filename = "DE output/C_d14C/DE_out.RData"
  filename = "DE output/C_d13C_d14C/DE_out.RData"
  load(file = filename)
  calibDummy <- outDE$optim$bestmem

  # Make sure the right calibration type is selected, based on which calibration was run
  # par@litterCalibration <- 0
  # par@soilCalibration <- 1
  
}

if(par@runSeries == 1){ # !!! MAKE SURE THE calibCheck PARAMETER IS CORRECT !!!!!!!
  # The file with parameter values is loaded
  load(file = "uniqueParam.RData") # This file contains the unique parameter values
                                   # for which the model needs to be run
  paramValues <- uniqueParam
  rm(uniqueParam)
  
  # The variable containing the outputs is created
  n <- nrow(paramValues) # The nuber of parameter combinations
  outputForSeries <- list("error" = rep(NA,n),
                          "Ctot" = matrix(NA, n, par@layerNumber),
                          "d13C" = matrix(NA, n, par@layerNumber),
                          "d14C" = matrix(NA, n, par@layerNumber),
                          "Cmic_R" = matrix(NA, n, par@layerNumber),
                          "Cmic_B" = matrix(NA, n, par@layerNumber),
                          "midDepth" = par@midDepth_layer,
                          "paramValues" = paramValues,
                          "relError_OC_POC" = rep(NA,n),
                          "relError_OC_MAOC" = rep(NA,n),
                          "relError_d13C_POC" = rep(NA,n),
                          "relError_d13C_MAOC" = rep(NA,n),
                          "relError_d14C_POC" = rep(NA,n),
                          "relError_d14C_MAOC" = rep(NA,n),
                          "POC_perc" = matrix(NA, n, par@layerNumber),
                          "bulkOC_perc" = matrix(NA, n, par@layerNumber),
                          "death_R_day" =  matrix(NA, n, par@layerNumber),
                          "death_B_day" =  matrix(NA, n, par@layerNumber))
  
  # Define which error measures need to be stored (based on which scenario is run)
  par@calibCheck <- data.frame(
    "C" = 1,
    "C_d13C" = 0,
    "C_d14C" = 0,
    "C_d13C_d14C" = 0)
  
}

# In case the parameter values for a SAFE sensitivity analysis are run
if(par@SAFE_sensitivity_isotopes == 1 | par@SAFE_sensitivity_allParameters == 1){
  
  # The parameter combinations, as defined by the Matlab version of SAFE
  paramValues <- read.csv("X_out.csv")

}

# ------------------------------------------------------------------------------
# Load the measured data
# ------------------------------------------------------------------------------

out <- loadMeasuredData(par@site)
measuredData_soil <- out$measuredData_soil
measuredData_litter <- out$measuredData_litter
rm(out)

# ------------------------------------------------------------------------------
# The input parameters are defined
# ------------------------------------------------------------------------------

out <- inputParameters(par) # Parameter values are defined in this function
par <- out$par
vegData_stage1 <- out$vegData_stage1 # The parameter values for vegetation stage 1
vegData_stage2 <- out$vegData_stage2 # The parameter values for vegetation stage 2
rm(out)

if(par@nVeg == 2){
  par@lastYearStage2 <- par@yearsStage2
}

# ------------------------------------------------------------------------------
# The model variables are initiated
# ------------------------------------------------------------------------------

out <- parameterInitialisation(par, vegData_stage1, vegData_stage2)
par = out$par                         # The model parameters
vegData_stage1 = out$vegData_stage1   # Vegetation data for stage 1
vegData_stage2 = out$vegData_stage2   # Vegetation data for stage 2
isotopicInputs = out$isotopicInputs   # Inputs for the different isotopes
isotopeRatios = out$isotopeRatios     # Ratios of the different isotopes
x_init = out$x_init                   # The initial size of the state variables
rm(out)

# ------------------------------------------------------------------------------
# The model is run
# ------------------------------------------------------------------------------

if(par@calibMode == 1){
  
  # -----------------------------
  # A parallel cluster is started
  # -----------------------------
  if(par@parallelCalib == 1){
    cl <- NCPUS(90)
    print(cl)
  }
  
  if(par@site == "Hainich"){
    
    if(par@litterCalibration == 1){
      
      # --------------------------------
      # Differential evolution algorithm
      # --------------------------------
      
      # The lower boundary
      lower <- c("KmB_POC_L_stage1" = 1e-2,
                 "KmB_DOC_L_stage1" = 1e-2,
                 "Vmax_L" = 1e-2)
      
      # The upper boundary
      upper <- c("KmB_POC_L_stage1" = 1e+3,
                 "KmB_DOC_L_stage1" = 1e+3,
                 "Vmax_L" = 1e2)
      
      # The variables that need to passed to DEoptim
      parVar <- list("par",
                     "vegData_stage1",
                     "vegData_stage2",
                     "isotopicInputs",
                     "isotopeRatios",
                     "x_init",
                     "measuredData_litter",
                     "measuredData_soil")
      
      outDE <- DEoptim(fn = runCarbonSimulation,
                       par1 = par,
                       vegData_stage1 = vegData_stage1,
                       vegData_stage2 = vegData_stage2,
                       isotopicInputs = isotopicInputs,
                       isotopeRatios = isotopeRatios,
                       measuredData_litter = measuredData_litter,
                       measuredData_soil = measuredData_soil,
                       x_init = x_init,
                       lower = lower,
                       upper = upper,
                       control = list(NP = 80,
                                      itermax = 30,
                                      cluster = cl,
                                      parVar = parVar)
      )
      
      save(outDE, file = "DE_out.RData")

    } else if (par@soilCalibration == 1){

      # --------------------------------
      # Differential evolution algorithm
      # --------------------------------
      
      # The lower boundary
      lower <- c("VmaxD_root_R_stage1" = 0.01,
                 "VmaxU_BioAv_stage1" = 0.01,
                 "Db0_stage1" = 1e-8,
                 "betaRoots" = 0.85,
                 "Km_ads" = 1e-8,
                 "Km_depol_M_stage1" = 1e-6,
                 "kDes_init" = 1e-6,
                 "VmaxD_M_stage1" = 0.1,
                 "Vmax_ads" = 0.1,
                 "Db_eFold_depth_stage1" = 0.01,
                 "advectionRate_polyC_stage1" = 0.01
                 )
      
      # The upper boundary
      upper <- c("VmaxD_root_R_stage1" = 1,
                 "VmaxU_BioAv_stage1" = 1,
                 "Db0_stage1" = 1e-4,
                 "betaRoots" = 0.97,
                 "Km_ads" = 1,
                 "Km_depol_M_stage1" = 1,
                 "kDes_init" = 1,
                 "VmaxD_M_stage1" = 1e3,
                 "Vmax_ads" = 1e3,
                 "Db_eFold_depth_stage1" = 0.5,
                 "advectionRate_polyC_stage1" = 1
                 )

      # The variables that need to passed to DEoptim
      parVar <- list("par",
                     "vegData_stage1",
                     "vegData_stage2",
                     "isotopicInputs",
                     "isotopeRatios",
                     "x_init",
                     "measuredData_litter",
                     "measuredData_soil")
      
      outDE <- DEoptim(fn = runCarbonSimulation,
                       par1 = par,
                       vegData_stage1 = vegData_stage1,
                       vegData_stage2 = vegData_stage2,
                       isotopicInputs = isotopicInputs,
                       isotopeRatios = isotopeRatios,
                       measuredData_litter = measuredData_litter,
                       measuredData_soil = measuredData_soil,
                       x_init = x_init,
                       lower = lower,
                       upper = upper,
                       control = list(NP = 180,
                                      itermax = 300,
                                      cluster = cl,
                                      parVar = parVar,
                                      storepopfrom = 1)
      )

      save(outDE, file = "DE_out.RData")
      
    }
    
  } # End if(site == "Hainich")
  
  # -------------------------------
  # The parallel cluster is stopped
  # -------------------------------
  if(par@parallelCalib == 1){
    stopCluster(cl)
  }

  
  # ----------------------------------------------------
  # If the model is run with the parameters for a SAFE sensitivity analysis for isotope parameters
  # ----------------------------------------------------
} else if (par@SAFE_sensitivity_isotopes == 1){

  # A cluster is started
  cl <- NCPUS(90)
  registerDoParallel(cl)
  print(cl)
  clusterExport(cl, c("par", "vegData_stage1", "vegData_stage2", "isotopicInputs",
                      "isotopeRatios", "x_init", "measuredData_litter", "measuredData_soil"))

  numRuns <- nrow(paramValues)

  output_SAFE <- matrix(NA,numRuns,6)

  parVar <- list("par",
                 "vegData_stage1",
                 "vegData_stage2",
                 "isotopicInputs",
                 "isotopeRatios",
                 "x_init",
                 "measuredData_litter",
                 "measuredData_soil")

  outForEach <- foreach(ii = 1:numRuns, .packages = c("deSolve", "psych")) %dopar% {

    output <- runCarbonSimulation(paramValues[ii,],par,vegData_stage1,vegData_stage2,
                               isotopicInputs,isotopeRatios,x_init,measuredData_litter,measuredData_soil)

    return(output)

  }

  stopCluster(cl)

  # The output from the foreach loop is reformatted
  for(jj in 1:numRuns){

    output_SAFE[jj,] <- outForEach[[jj]]

  }

  save(outForEach, file = "outForEach.RData")
  save(output_SAFE, file = "output_SAFE.RData")
  
    
  # ----------------------------------------------------
  # If the model is run with the parameters for a SAFE sensitivity analysis for all parameters
  # ----------------------------------------------------
} else if (par@SAFE_sensitivity_allParameters == 1){
  
  # A cluster is started
  cl <- NCPUS(90)
  registerDoParallel(cl)
  print(cl)
  clusterExport(cl, c("par", "vegData_stage1", "vegData_stage2", "isotopicInputs",
                      "isotopeRatios", "x_init", "measuredData_litter", "measuredData_soil"))
  
  numRuns <- nrow(paramValues)
  
  output_SAFE <- matrix(NA,numRuns,7)
  
  parVar <- list("par",
                 "vegData_stage1",
                 "vegData_stage2",
                 "isotopicInputs",
                 "isotopeRatios",
                 "x_init",
                 "measuredData_litter",
                 "measuredData_soil")
  
  outForEach <- foreach(ii = 1:numRuns, .packages = c("deSolve", "psych")) %dopar% {
    
    output <- runCarbonSimulation(paramValues[ii,],par,vegData_stage1,vegData_stage2,
                                  isotopicInputs,isotopeRatios,x_init,measuredData_litter,measuredData_soil)
    
    return(output)
    
  }
  
  stopCluster(cl)
  
  # The output from the foreach loop is reformatted
  for(jj in 1:numRuns){
    
    output_SAFE[jj,] <- outForEach[[jj]]
    
  }
  
  save(outForEach, file = "outForEach.RData")
  save(output_SAFE, file = "output_SAFE.RData")
  
  # ----------------------------------------------------
  # If the model is run for a series of parameter values
  # ----------------------------------------------------
  
}else if (par@runSeries == 1){

  # A cluster is started
  cl <- NCPUS(90)
  registerDoParallel(cl)
  print(cl)
  clusterExport(cl, c("par", "vegData_stage1", "vegData_stage2", "isotopicInputs",
                                          "isotopeRatios", "x_init", "measuredData_litter", "measuredData_soil"))

  numRuns <- nrow(paramValues)

  parVar <- list("par",
                 "vegData_stage1",
                 "vegData_stage2",
                 "isotopicInputs",
                 "isotopeRatios",
                 "x_init",
                 "measuredData_litter",
                 "measuredData_soil")

  outForEach <- foreach(ii = 1:numRuns, .packages = c("deSolve", "psych")) %dopar% {

  out <- runCarbonSimulation(paramValues[ii,],par,vegData_stage1,vegData_stage2,
                            isotopicInputs,isotopeRatios,x_init,measuredData_litter,measuredData_soil)
  
  endOf <- nrow(out$res$C_soil)
  
  out_C <- (out$res$C_soil[endOf,] / par@soilMassArray) * 100
  out_d13C <- out$d13C$soilC[endOf,]
  out_d14C <- out$d14C$soilC[endOf,]
  out_Cmic_R <- out$res$MIC_R[endOf,]
  out_Cmic_B <- out$res$MIC_B[endOf,]
  
  error <- out$error
  relError_OC_POC <- out$relError_OC_POC
  relError_OC_MAOC <- out$relError_OC_MAOC
  relError_d13C_POC <- out$relError_d13C_POC
  relError_d13C_MAOC <- out$relError_d13C_MAOC
  relError_d14C_POC <- out$relError_d14C_POC
  relError_d14C_MAOC <- out$relError_d14C_MAOC
  
  POC_perc <- (out$res$POC_R[endOf,] / par@soilMassArray) * 100
  bulkOC_perc <- ((out$res$POLY_B[endOf,] + out$res$MIN_B[endOf,] + out$res$MIC_B[endOf,]) / par@soilMassArray) * 100
  
  death_R_day <- out$res$death_R_day
  death_B_day <- out$res$death_B_day

  return(list("C" = out_C,
              "d13C" = out_d13C,
              "d14C" = out_d14C,
              "Cmic_R" = out_Cmic_R,
              "Cmic_B" = out_Cmic_B,
              "error" = error,
              "relError_OC_POC" = relError_OC_POC,
              "relError_OC_MAOC" = relError_OC_MAOC,
              "relError_d13C_POC" = relError_d13C_POC,
              "relError_d13C_MAOC" = relError_d13C_MAOC,
              "relError_d14C_POC" = relError_d14C_POC,
              "relError_d14C_MAOC" = relError_d14C_MAOC,
              "POC_perc" = POC_perc,
              "bulkOC_perc" = bulkOC_perc,
              "death_R_day" = death_R_day,
              "death_B_day" = death_B_day))
  }
  
  stopCluster(cl)
  
  # The output from the foreach loop is reformatted
  for(jj in 1:numRuns){
    
    outputForSeries$Ctot[jj,] <- outForEach[[jj]]$C
    outputForSeries$d13C[jj,] <- outForEach[[jj]]$d13C
    outputForSeries$d14C[jj,] <- outForEach[[jj]]$d14C
    outputForSeries$Cmic_R[jj,] <- outForEach[[jj]]$Cmic_R
    outputForSeries$Cmic_B[jj,] <- outForEach[[jj]]$Cmic_B
    outputForSeries$error[jj] <- outForEach[[jj]]$error
    outputForSeries$relError_OC_POC[jj] <- outForEach[[jj]]$relError_OC_POC
    outputForSeries$relError_OC_MAOC[jj] <- outForEach[[jj]]$relError_OC_MAOC
    outputForSeries$relError_d13C_POC[jj] <- outForEach[[jj]]$relError_d13C_POC
    outputForSeries$relError_d13C_MAOC[jj] <- outForEach[[jj]]$relError_d13C_MAOC
    outputForSeries$relError_d14C_POC[jj] <- outForEach[[jj]]$relError_d14C_POC
    outputForSeries$relError_d14C_MAOC[jj] <- outForEach[[jj]]$relError_d14C_MAOC
    outputForSeries$POC_perc[jj,] <- outForEach[[jj]]$POC_perc
    outputForSeries$bulkOC_perc[jj,] <- outForEach[[jj]]$bulkOC_perc
    
  }
  
  # The output is saved
  if(par@calibCheck$C == 1){
    save(outputForSeries, file = "outputForSeries_C.RData")
  } else if(par@calibCheck$C_d13C == 1){
    save(outputForSeries, file = "outputForSeries_C_d13C.RData")
  } else if(par@calibCheck$C_d14C == 1){
    save(outputForSeries, file = "outputForSeries_C_d14C.RData")
  } else if(par@calibCheck$C_d13C_d14C == 1){
    save(outputForSeries, file = "outputForSeries_c_d13C_d14C.RData")
  }

  # ----------------------------------------------------
  # If a single model run is performed
  # ----------------------------------------------------
  
} else if (par@normalRun == 1) {

  # A dummy variable for the calibration parameters
  if(par@useCalibrationResults == 0){
    calibDummy = NA
  }
  
  out <- runCarbonSimulation(calibDummy,par, vegData_stage1,vegData_stage2,
                             isotopicInputs,isotopeRatios,x_init,measuredData_litter,measuredData_soil)
  
  # The output is formatted
  res <- out$res
  d13C <- out$d13C
  d14C <- out$d14C
  rm(out)
   
}
  
# ------------------------------------------------------------------------------
# The results are plotted
# ------------------------------------------------------------------------------

if(par@calibMode == 0 & par@runSeries == 0){
  
  out <- plotting(res, d13C, d14C, par)
  
  pl1 <- out[[1]]
  pl2 <- out[[2]]
  pl3 <- out[[3]]
  
  # quartz(title = "litter carbon", width = 8, height = 5)
  # ggarrange(pl1, pl2, pl3, nrow = 3, labels = c("(A)", "(B)", "(C)"), font.label=list(color="black",size=12))
  
  pc1 <- out[[4]]
  pc2 <- out[[5]]
  pc3 <- out[[6]]
  
  quartz(title = "Depth profiles", width = 16*0.39, height = 7*0.39)
  ggarrange(pc1,pc2,pc3, nrow = 1, labels = c("(A)", "(B)", "(C)"), font.label=list(color="black",size=12))
  
  # Turnover time of microbes
  p_TT <- out[[7]]
  quartz(title = "Turnover time", width = 16*0.39, height = 16*0.39)
  p_TT
  
}

# ---------------------------------
# The time the run took is recorded
# ---------------------------------

endTime <- proc.time()
print(elapsedTime <-endTime[3] - startTime[3])